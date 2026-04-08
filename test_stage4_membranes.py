"""Stage 4 verification: membrane self-assembly + compartment finder.

Three checks:
  1. Regression: membrane_threshold_high=999 (default) → no membrane (Stage 3 byte-compatible).
  2. Synthetic: hand-place a 5×5 lipid blob → surface-only rule produces a hollow ring,
     and find_compartments finds the 3×3 interior.
  3. Emergent: hand-seed a strong lipid catalyst → run the simulation → eventually
     a closed compartment forms. Sweep lipid_chemotaxis to find the best regime.
"""
import numpy as np
from sim import World, WorldConfig, Strand
from analysis.compartments import find_compartments, compartment_stats


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
    )
    base.update(overrides)
    return WorldConfig(**base)


def main():
    print("─── Test 1: REGRESSION (membrane_threshold_high=999) ───")
    cfg = make_cfg(p_mut=0.01, k_lipid=0.5)
    w = World(cfg, seed=42)
    # Seed a lipid catalyst so lipids actually accumulate
    seed = Strand([0, 3, 0, 1, 2], x=12, y=12, gen=0)
    w.strands.append(seed)
    for t in range(300):
        w.step(env_E=0.25)
    n_membrane = int(w.is_membrane.sum())
    print(f"  lipid total: {int(w.lipid.sum())}  is_membrane cells: {n_membrane}")
    assert n_membrane == 0, "REGRESSION FAILED: membrane appeared with default threshold"
    print("  ✓ no membranes when threshold disabled\n")

    print("─── Test 2: SYNTHETIC blob → ring topology ───")
    cfg = make_cfg(membrane_threshold_high=8, membrane_threshold_low=4)
    w = World(cfg, seed=0)
    # Hand-place a 5×5 solid lipid blob centered at (10,10)
    for i in range(8, 13):
        for j in range(8, 13):
            w.lipid[i, j] = 20  # well above HIGH=8
    # Run one step to update membrane state (no other dynamics needed)
    w.step(env_E=0.25)
    print(f"  lipid total: {int(w.lipid.sum())}  is_membrane cells: {int(w.is_membrane.sum())}")
    # Pretty-print the membrane region around the blob
    print("  is_membrane (M=membrane, I=lipid_filled-but-not-membrane, .=empty):")
    for i in range(6, 15):
        row = "    "
        for j in range(6, 15):
            if w.is_membrane[i, j]:
                row += "M "
            elif w._lipid_filled[i, j]:
                row += "I "
            else:
                row += ". "
        print(row)
    exterior, interiors = find_compartments(w.is_membrane)
    stats = compartment_stats(exterior, interiors)
    print(f"  compartments: {stats['n_compartments']}  largest_interior: {stats['largest_interior']}")
    assert stats["n_compartments"] >= 1, "SYNTHETIC FAILED: no compartment from solid blob"
    assert stats["largest_interior"] >= 9, (
        f"SYNTHETIC FAILED: expected interior of ≥9 (3×3), got {stats['largest_interior']}")
    print("  ✓ surface-only rule produces ring + interior\n")

    print("─── Test 3: EMERGENT compartment formation (sweep chemotaxis) ───")
    print("  Hand-seed a strong lipid catalyst, run 1500 steps, check compartments.")
    print(f"  {'chemotaxis':>11} {'max_compartments':>17} {'final_n_compartments':>22} {'max_lipid':>10}")
    best_chemo = None
    best_max = 0
    for chemo in [0.0, 0.3, 0.5, 0.7, 0.9]:
        cfg = make_cfg(
            spont_rate=0.0, extend_rate=0.0,  # isolate the seeded catalyst
            k_lipid=1.0,                       # strong synthesis
            lipid_chemotaxis=chemo,
            membrane_threshold_high=8,
            membrane_threshold_low=4,
        )
        w = World(cfg, seed=42)
        # Seed several catalyst strands at different cells so multiple lipid sources
        # exist (encourages multi-blob → multi-compartment dynamics)
        for (x, y) in [(6, 6), (6, 18), (18, 6), (18, 18), (12, 12)]:
            seed = Strand([0, 3, 0, 0, 3, 0], x=x, y=y, gen=0)  # double motif!
            w.strands.append(seed)

        max_compartments = 0
        max_lipid = 0
        for t in range(1500):
            w.step(env_E=0.25)
            if (t + 1) % 50 == 0:
                _, interiors = find_compartments(w.is_membrane)
                if len(interiors) > max_compartments:
                    max_compartments = len(interiors)
                if int(w.lipid.max()) > max_lipid:
                    max_lipid = int(w.lipid.max())
        _, final_interiors = find_compartments(w.is_membrane)
        print(f"  {chemo:>11.1f} {max_compartments:>17d} {len(final_interiors):>22d} {max_lipid:>10d}")
        if max_compartments > best_max:
            best_max = max_compartments
            best_chemo = chemo

    print(f"\n  Best chemotaxis: {best_chemo}  (max simultaneous compartments: {best_max})")
    if best_max >= 1:
        print("  ✓ at least one emergent compartment formed in some regime")
    else:
        print("  ⚠ no emergent compartment formed — physics tuning needed")

    print("\nStage 4 done.")


if __name__ == "__main__":
    main()
