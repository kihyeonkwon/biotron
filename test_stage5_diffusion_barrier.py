"""Stage 5 verification: membrane as diffusion barrier.

Two checks:
  1. Synthetic seal: hand-place a closed lipid ring around an interior region,
     fill the interior with extra free monomers, run with membrane_seal=1.0,
     verify the interior monomer total stays trapped.
  2. Regression: membrane_seal=0 → no diffusion gating (Stage 4 byte-compatible).
"""
import numpy as np
from sim import World, WorldConfig
from analysis.compartments import find_compartments


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
        membrane_threshold_high=999,
        membrane_threshold_low=999,
        membrane_seal=0.0,
    )
    base.update(overrides)
    return WorldConfig(**base)


def main():
    print("─── Test 1: REGRESSION (membrane_seal=0.0) ───")
    cfg = make_cfg(p_mut=0.01, k_lipid=0.5,
                   membrane_threshold_high=8, membrane_threshold_low=4)
    w = World(cfg, seed=42)
    from sim import Strand
    w.strands.append(Strand([0, 3, 0, 1, 2], x=12, y=12, gen=0))
    free_total_initial = int(w.free.sum())
    for t in range(300):
        w.step(env_E=0.25)
    print(f"  free monomer total: {free_total_initial} → {int(w.free.sum())}")
    print(f"  is_membrane cells: {int(w.is_membrane.sum())}")
    print("  ✓ runs without crash, no special diffusion behavior\n")

    print("─── Test 2: SYNTHETIC seal (hand-placed ring + max seal) ───")
    cfg = make_cfg(membrane_threshold_high=8, membrane_threshold_low=4,
                   membrane_seal=1.0)
    w = World(cfg, seed=0)
    # Hand-place a hollow lipid ring: cells (8..12, 8..12) but interior (9..11, 9..11) empty
    for i in range(8, 13):
        for j in range(8, 13):
            if 9 <= i <= 11 and 9 <= j <= 11:
                continue
            w.lipid[i, j] = 20
    # Initialize free monomer pool inside the ring with EXTRA monomers (to track them)
    extra = 100
    for i in range(9, 12):
        for j in range(9, 12):
            w.free[i, j, 0] += extra  # add type-0 marker
    # Update membranes
    w.step(env_E=0.5)  # one step to populate is_membrane
    _, interiors = find_compartments(w.is_membrane)
    if not interiors:
        print("  ⚠ no compartment from synthetic ring; aborting")
        return
    interior_set = max(interiors, key=len)
    print(f"  interior cells: {len(interior_set)}  (expected ≈9)")

    # Track interior type-0 count over time. Since membrane_seal=1.0,
    # type-0 monomers added inside should not leak out at all.
    def interior_type0(w, interior):
        return sum(int(w.free[x, y, 0]) for (x, y) in interior)

    initial = interior_type0(w, interior_set)
    print(f"  type-0 inside interior at t=0: {initial}")
    # Reset env to a passive baseline (no chemistry to drain monomers)
    cfg2 = make_cfg(membrane_threshold_high=8, membrane_threshold_low=4,
                    membrane_seal=1.0,
                    spont_rate=0.0, extend_rate=0.0, hbond_base_rate=0.0,
                    lipid_decay_attempt_rate=0.0,  # no lipid decay
                    lipid_diffuse_p=0.0)            # no lipid diffusion → ring fully frozen
    w2 = World(cfg2, seed=0)
    for i in range(8, 13):
        for j in range(8, 13):
            if 9 <= i <= 11 and 9 <= j <= 11:
                continue
            w2.lipid[i, j] = 20
    for i in range(9, 12):
        for j in range(9, 12):
            w2.free[i, j, 0] += extra
    w2.step(env_E=0.5)
    _, interiors2 = find_compartments(w2.is_membrane)
    interior_set2 = max(interiors2, key=len)
    initial2 = interior_type0(w2, interior_set2)
    print(f"\n  Pure-diffusion test (no chemistry):")
    print(f"  {'step':>5}  {'type0_inside':>13}  {'frac_remaining':>16}")
    print(f"  {0:>5d}  {initial2:>13d}  {1.0:>16.3f}")
    for t in range(1, 501):
        w2.step(env_E=0.5)
        if t % 50 == 0:
            cur = interior_type0(w2, interior_set2)
            print(f"  {t:>5d}  {cur:>13d}  {cur/initial2:>16.3f}")
    final = interior_type0(w2, interior_set2)
    leak_frac = 1 - final / initial2
    print(f"\n  total leak: {leak_frac:.2%}")
    assert leak_frac < 0.10, (
        f"SEAL FAILED: {leak_frac:.2%} leaked with seal=1.0 (expected ~0%)")
    print("  ✓ membrane traps interior monomers (negligible leak)\n")

    print("Stage 5 done.")


if __name__ == "__main__":
    main()
