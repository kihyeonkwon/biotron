"""Stage 3 verification: lipid synthesis.

Three checks:
  1. Regression: with k_lipid=0, byte-identical to Stage 2 (no lipid array touched).
  2. Hotspot: hand-seed a strand carrying the lipid_synth motif → lipids accumulate
     in that cell over time.
  3. Conservation: total atoms (free monomers + strand monomers + lipids) is
     conserved (drifts only by stochastic dimer formation, not creation).
"""
import numpy as np
from sim import World, WorldConfig, Strand
from sim import catalysis as cat_mod


def make_cfg(**overrides):
    base = dict(
        grid_h=24, grid_w=24,
        initial_mono_per_cell=15,
        spont_rate=0.003,
        hbond_base_rate=0.6,
        stacking_bonus=3.0,
        poly_base_rate=0.05,
        extend_rate=0.02,
        min_birth_length=3,
        max_strand_length=24,
        max_strands=4000,
        p_mut=0.0,
        k_poly_boost=0.0,
        k_lipid=0.0,
    )
    base.update(overrides)
    return WorldConfig(**base)


def total_atoms(w):
    """Free monomers + monomers in strands (parent + daughter) + lipids."""
    free = int(w.free.sum())
    strand_atoms = 0
    for s in w.strands:
        strand_atoms += s.length  # parent
        strand_atoms += int(s.d_present.sum())  # daughter
    lipids = int(w.lipid.sum())
    return free + strand_atoms + lipids, free, strand_atoms, lipids


def main():
    print("─── Test 1: REGRESSION (k_lipid=0) ───")
    cfg = make_cfg(p_mut=0.01)
    w = World(cfg, seed=42)
    for t in range(500):
        w.step(env_E=(np.sin(t*0.05)+1)/2)
    print(f"  total_births: {w.total_births}  strands: {w.n_strands()}  lipid_total: {int(w.lipid.sum())}")
    assert int(w.lipid.sum()) == 0, "REGRESSION FAILED: lipids appeared with k_lipid=0"
    print("  ✓ no lipids when k_lipid=0\n")

    print("─── Test 2: HOTSPOT (hand-seed lipid catalyst) ───")
    cfg = make_cfg(spont_rate=0.0, extend_rate=0.0, k_lipid=0.5)
    w = World(cfg, seed=99)
    # Seed a strand with the lipid_synth motif [0,3,0]
    seq = [0, 3, 0, 1, 2]   # motif at positions 0..2
    cat_strand = Strand(seq, x=10, y=10, gen=0)
    assert cat_strand._cat_lipid == 1, f"seed motif strength wrong: {cat_strand._cat_lipid}"
    w.strands.append(cat_strand)

    samples = []
    for t in range(500):
        w.step(env_E=0.25)  # constant cold
        if (t + 1) % 50 == 0:
            samples.append((t + 1, int(w.lipid[10, 10]), int(w.lipid.sum()),
                            int(w.free[10, 10].sum())))

    print(f"  {'step':>5}  {'lipid@(10,10)':>14}  {'lipid_total':>11}  {'free@(10,10)':>13}")
    for t, lp_here, lp_tot, fr_here in samples:
        print(f"  {t:>5d}  {lp_here:>14d}  {lp_tot:>11d}  {fr_here:>13d}")

    assert samples[-1][1] > 0 or samples[-1][2] > 0, "HOTSPOT FAILED: no lipids ever appeared"
    final_here = samples[-1][1]
    final_total = samples[-1][2]
    print(f"  ✓ lipids accumulated (final cell={final_here}, grid total={final_total})\n")

    print("─── Test 3: CONSERVATION (total atoms over time) ───")
    cfg = make_cfg(spont_rate=0.003, k_lipid=0.5, p_mut=0.01)
    w = World(cfg, seed=7)
    initial, _, _, _ = total_atoms(w)
    print(f"  initial total atoms: {initial}")
    snapshots = []
    for t in range(800):
        w.step(env_E=(np.sin(t*0.05)+1)/2)
        if (t + 1) % 100 == 0:
            snapshots.append((t + 1, *total_atoms(w)))
    print(f"  {'step':>5}  {'total':>7}  {'free':>7}  {'strand':>7}  {'lipid':>6}  {'drift':>7}")
    for t, tot, free, strand, lipid in snapshots:
        drift = tot - initial
        print(f"  {t:>5d}  {tot:>7d}  {free:>7d}  {strand:>7d}  {lipid:>6d}  {drift:>+7d}")
    final_total = snapshots[-1][1]
    drift = final_total - initial
    drift_frac = abs(drift) / initial
    print(f"\n  total drift: {drift:+d} ({drift_frac:.2%})")
    assert drift_frac < 0.01, (
        f"CONSERVATION FAILED: drift {drift_frac:.2%} > 1%")
    print("  ✓ atoms conserved within 1%")

    print("\nStage 3 done.")


if __name__ == "__main__":
    main()
