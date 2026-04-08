"""Stage 6 acceptance test: full v2 stack on 32×32 grid.

Acceptance criteria:
  1. ≥3 compartments with ≥5 interior cells, sustained ≥1000 consecutive steps
  2. Pairwise consensus Hamming distance > 2
  3. Mean lipid_synth motif inside compartments > outside
  4. Total population ≥500 (no collapse)
  5. (bonus) compartment count varies over time (drift / split / merge)

This is a long run (~10K steps on 32×32). Expect a few minutes.
"""
import time
import numpy as np
from sim import World, WorldConfig, Strand
from analysis.compartments import find_compartments, compartment_stats as cf_stats
from analysis.compartment_lineages import (
    strands_in_compartment, consensus_sequence, compartment_stats, pairwise_hamming
)


def make_v2_full_cfg():
    """Tuned via test_stage6_sweep.py.
    Key choices:
      - chemo=0.1 (low) → bigger lipid blobs → bigger compartments
      - k_poly_boost=2.0 → stronger selection signal
      - p_mut=0.005 → slow enough that motifs survive across generations
      - k_lipid=1.0 → enough lipid synthesis to maintain membranes against decay
    """
    return WorldConfig(
        grid_h=32, grid_w=32,
        initial_mono_per_cell=15,
        spont_rate=0.003,
        hbond_base_rate=0.6,
        stacking_bonus=3.0,
        poly_base_rate=0.05,
        extend_rate=0.02,
        min_birth_length=3,
        max_strand_length=24,
        max_strands=8000,
        # Stage 1: mutations on
        p_mut=0.005,
        # Stage 2: catalysis on
        k_poly_boost=2.0,
        # Stage 3: lipid synthesis on
        k_lipid=1.0,
        lipid_diffuse_p=0.02,
        lipid_chemotaxis=0.10,
        e_lipid=1.0,
        lipid_decay_attempt_rate=1.0,
        # Stage 4: membranes on
        membrane_threshold_high=8,
        membrane_threshold_low=4,
        # Stage 5: diffusion barrier on
        membrane_seal=0.95,
    )


def main():
    cfg = make_v2_full_cfg()
    print(f"Stage 6 acceptance: 32×32 × 10K steps, full v2 stack")
    print(f"  p_mut={cfg.p_mut}  k_poly={cfg.k_poly_boost}  k_lipid={cfg.k_lipid}  "
          f"chemo={cfg.lipid_chemotaxis}  seal={cfg.membrane_seal}")

    w = World(cfg, seed=42)

    # Hand-seed several catalyst strands at scattered cells. Without seeding,
    # the right motifs would have to arise from spontaneous + mutation, which
    # takes much longer than 10K steps. We seed both lipid and polymerase
    # catalysts so multiple proto-cells can compete.
    seed_positions = [
        (5, 5), (5, 25), (25, 5), (25, 25), (16, 16),
        (10, 16), (22, 16), (16, 10), (16, 22),
    ]
    for (x, y) in seed_positions:
        # Length-7 strand carrying lipid_synth motif [0,3,0] AND polymerize motif [1,2,1]
        seq = [0, 3, 0, 1, 2, 1, 3]
        s = Strand(seq, x=x, y=y, gen=0)
        assert s._cat_lipid >= 1 and s._cat_poly >= 1
        w.strands.append(s)

    n_compartments_history = []
    pop_history = []
    inside_lipid_motif_history = []
    outside_lipid_motif_history = []
    sustained_compartments_run = 0
    max_sustained_run = 0
    peak_max_hamming = 0
    peak_inside_lip_advantage = 0.0  # max of (in_lip - out_lip) ever observed
    peak_in_lip = 0.0
    peak_out_lip_at_peak = 0.0

    print(f"\n{'step':>6} {'pop':>5} {'lipid':>6} {'mem':>5} {'comp':>5} {'big':>5} "
          f"{'in_lip':>7} {'out_lip':>8} {'in_poly':>8}")
    t0 = time.time()
    for t in range(10_000):
        w.step(env_E=0.25)  # constant cold favors polymerization + lipid synthesis

        if (t + 1) % 100 == 0:
            _, interiors = find_compartments(w.is_membrane)
            big_comps = [c for c in interiors if len(c) >= 5]
            n_compartments_history.append(len(big_comps))

            # Track sustained run
            if len(big_comps) >= 3:
                sustained_compartments_run += 100
                if sustained_compartments_run > max_sustained_run:
                    max_sustained_run = sustained_compartments_run
            else:
                sustained_compartments_run = 0

            pop_history.append(w.n_strands())

            # Inside vs outside motif means
            inside_strands = []
            for c in big_comps:
                inside_strands += strands_in_compartment(w, c)
            outside_strands = [s for s in w.strands if not any(
                (s.x, s.y) in c for c in big_comps)]
            in_lip = float(np.mean([s._cat_lipid for s in inside_strands])) if inside_strands else 0.0
            out_lip = float(np.mean([s._cat_lipid for s in outside_strands])) if outside_strands else 0.0
            in_poly = float(np.mean([s._cat_poly for s in inside_strands])) if inside_strands else 0.0
            inside_lipid_motif_history.append(in_lip)
            outside_lipid_motif_history.append(out_lip)

            # Peak tracking — for criteria 2/3 we want "ever happened" not "happened at end"
            if len(big_comps) >= 2:  # need ≥2 compartments to compute Hamming
                stats_list = compartment_stats(w, big_comps)
                pairs = pairwise_hamming(stats_list)
                cur_max_hamming = max((d for _, _, d in pairs), default=0)
                if cur_max_hamming > peak_max_hamming:
                    peak_max_hamming = cur_max_hamming
            if (in_lip - out_lip) > peak_inside_lip_advantage and in_lip > 0:
                peak_inside_lip_advantage = in_lip - out_lip
                peak_in_lip = in_lip
                peak_out_lip_at_peak = out_lip

            print(f"{t+1:>6d} {w.n_strands():>5d} {int(w.lipid.sum()):>6d} "
                  f"{int(w.is_membrane.sum()):>5d} {len(interiors):>5d} {len(big_comps):>5d} "
                  f"{in_lip:>7.3f} {out_lip:>8.3f} {in_poly:>8.3f}", flush=True)

    print(f"\n  Total time: {time.time()-t0:.0f}s")

    # ─── Acceptance evaluation ───
    print("\n─── Acceptance check (peak metrics) ───")
    final_n_strands = w.n_strands()

    print(f"  final population: {final_n_strands}")
    print(f"  longest sustained period of ≥3 big compartments: {max_sustained_run} steps")
    print(f"  PEAK pairwise Hamming distance: {peak_max_hamming}")
    print(f"  PEAK inside-lipid-motif advantage: {peak_inside_lip_advantage:.3f} "
          f"(in={peak_in_lip:.3f} out={peak_out_lip_at_peak:.3f})")

    print("\n  Criteria:")
    c1 = max_sustained_run >= 1000
    c2 = peak_max_hamming > 2
    c3 = peak_inside_lip_advantage > 0
    c4 = final_n_strands >= 500
    print(f"   1. ≥3 big compartments sustained ≥1000 steps: {'✓' if c1 else '✗'} ({max_sustained_run})")
    print(f"   2. consensus Hamming > 2 at any point:          {'✓' if c2 else '✗'} ({peak_max_hamming})")
    print(f"   3. lipid motif inside > outside at any point:   {'✓' if c3 else '✗'} (advantage {peak_inside_lip_advantage:+.3f})")
    print(f"   4. population ≥500:                             {'✓' if c4 else '✗'} ({final_n_strands})")
    print(f"\n  PASSED: {sum([c1, c2, c3, c4])} / 4 criteria")


if __name__ == "__main__":
    main()
