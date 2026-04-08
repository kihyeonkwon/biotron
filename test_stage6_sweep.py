"""Quick parameter sweep for Stage 6 acceptance criteria.

Tests several combinations of (k_poly_boost, p_mut, lipid_chemotaxis, k_lipid)
on shorter (3000-step) runs to find a regime that hits more criteria.
"""
import time
import numpy as np
from sim import World, WorldConfig, Strand
from analysis.compartments import find_compartments
from analysis.compartment_lineages import strands_in_compartment, compartment_stats, pairwise_hamming


def make_cfg(k_poly, p_mut, chemo, k_lipid):
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
        p_mut=p_mut,
        k_poly_boost=k_poly,
        k_lipid=k_lipid,
        lipid_diffuse_p=0.02,
        lipid_chemotaxis=chemo,
        e_lipid=1.0,
        lipid_decay_attempt_rate=1.0,
        membrane_threshold_high=8,
        membrane_threshold_low=4,
        membrane_seal=0.95,
    )


def run_one(k_poly, p_mut, chemo, k_lipid, n_steps=3000):
    cfg = make_cfg(k_poly, p_mut, chemo, k_lipid)
    w = World(cfg, seed=42)
    seed_seqs = [
        [0, 3, 0, 1, 2, 1, 3],   # both motifs
        [0, 3, 0, 0, 3, 0, 1],   # double lipid motif
        [1, 2, 1, 0, 3, 0, 2],   # both
        [0, 3, 0, 1, 2, 1, 3],   # both
    ]
    seed_positions = [
        (5, 5), (5, 25), (25, 5), (25, 25),
        (16, 16), (10, 16), (22, 16), (16, 10), (16, 22),
        (8, 8), (8, 24), (24, 8), (24, 24),
    ]
    for i, (x, y) in enumerate(seed_positions):
        seq = seed_seqs[i % len(seed_seqs)]
        w.strands.append(Strand(seq, x=x, y=y, gen=0))

    max_big = 0
    sustained_big_run = 0
    max_sustained = 0
    final_in_lip = 0
    final_out_lip = 0
    final_n_compartments = 0
    final_max_hamming = 0
    for t in range(n_steps):
        w.step(env_E=0.25)
        if (t + 1) % 100 == 0:
            _, interiors = find_compartments(w.is_membrane)
            big = [c for c in interiors if len(c) >= 5]
            if len(big) > max_big:
                max_big = len(big)
            if len(big) >= 3:
                sustained_big_run += 100
                if sustained_big_run > max_sustained:
                    max_sustained = sustained_big_run
            else:
                sustained_big_run = 0
    # Final
    _, interiors = find_compartments(w.is_membrane)
    big = [c for c in interiors if len(c) >= 5]
    final_n_compartments = len(big)
    inside_strands = []
    for c in big:
        inside_strands += strands_in_compartment(w, c)
    outside_strands = [s for s in w.strands if not any((s.x, s.y) in c for c in big)]
    final_in_lip = float(np.mean([s._cat_lipid for s in inside_strands])) if inside_strands else 0.0
    final_out_lip = float(np.mean([s._cat_lipid for s in outside_strands])) if outside_strands else 0.0
    stats_list = compartment_stats(w, big)
    pairs = pairwise_hamming(stats_list)
    final_max_hamming = max((d for _, _, d in pairs), default=0)
    return {
        "max_big": max_big,
        "max_sustained": max_sustained,
        "final_big": final_n_compartments,
        "final_in_lip": final_in_lip,
        "final_out_lip": final_out_lip,
        "final_pop": w.n_strands(),
        "final_max_hamming": final_max_hamming,
        "final_lipid": int(w.lipid.sum()),
    }


def main():
    configs = [
        # name, k_poly, p_mut, chemo, k_lipid
        ("baseline",     0.5, 0.010, 0.30, 0.5),
        ("strong_sel",   2.0, 0.005, 0.30, 0.5),
        ("strong_sel2",  5.0, 0.005, 0.30, 1.0),
        ("low_chemo",    2.0, 0.005, 0.10, 1.0),
        ("aggressive",   5.0, 0.002, 0.10, 1.5),
        ("low_mut",      2.0, 0.001, 0.30, 1.0),
    ]
    print(f"{'name':>14} {'maxBig':>8} {'sustain':>9} {'finBig':>8} "
          f"{'inLip':>7} {'outLip':>8} {'maxHam':>8} {'pop':>5} {'lipid':>6} {'time':>6}")
    print("─" * 95)
    for name, k_poly, p_mut, chemo, k_lipid in configs:
        t0 = time.time()
        r = run_one(k_poly, p_mut, chemo, k_lipid, n_steps=3000)
        dt = time.time() - t0
        print(f"{name:>14} {r['max_big']:>8d} {r['max_sustained']:>9d} "
              f"{r['final_big']:>8d} {r['final_in_lip']:>7.2f} "
              f"{r['final_out_lip']:>8.3f} {r['final_max_hamming']:>8d} "
              f"{r['final_pop']:>5d} {r['final_lipid']:>6d} {dt:>5.0f}s", flush=True)


if __name__ == "__main__":
    main()
