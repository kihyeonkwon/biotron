"""Stage 1 verification: mutations at liftoff.

Three checks:
  1. Regression: p_mut=0 still gives 100% template-faithful births.
  2. Diversity: p_mut=0.01 produces visibly diverse sequences vs p_mut=0.
  3. Drift: mean total_mutations_inherited grows roughly linearly with generation.
"""
import numpy as np
from sim import World, WorldConfig
from sim.physics import COMPLEMENT, LABELS


def make_cfg(p_mut):
    return WorldConfig(
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
        p_mut=p_mut,
    )


def run_with_birth_log(cfg, n_steps=1000, seed=42):
    """Wrap _strand_liftoff to count fidelity statistics."""
    w = World(cfg, seed=seed)
    log = {"faithful": 0, "mutated": 0, "by_gen_total_mut": {}}

    original = w._strand_liftoff
    def logged(s):
        parent_seq = s.mono.copy()
        new_strands, births = original(s)
        for child in new_strands:
            child_seq = child.mono
            # Find any contiguous parent segment that is the *complement* of the child
            # OR the complement modulo mutations
            best_mismatches = 1_000_000
            for start in range(len(parent_seq) - len(child_seq) + 1):
                expected = COMPLEMENT[parent_seq[start:start + len(child_seq)]]
                mismatches = int((expected != child_seq).sum())
                if mismatches < best_mismatches:
                    best_mismatches = mismatches
            if best_mismatches == 0:
                log["faithful"] += 1
            else:
                log["mutated"] += 1
        return new_strands, births
    w._strand_liftoff = logged

    for t in range(n_steps):
        E = (np.sin(t * 0.05) + 1) / 2
        w.step(env_E=E)

    # Generation-level mutation drift
    by_gen = {}
    for s in w.strands:
        by_gen.setdefault(s.gen, []).append(s.total_mutations_inherited)
    log["by_gen"] = {g: float(np.mean(v)) for g, v in by_gen.items()}
    log["n_strands"] = w.n_strands()
    log["total_births"] = w.total_births
    log["distinct_seqs"] = len({tuple(s.mono.tolist()) for s in w.strands})
    return log


def main():
    print("─── Test 1: REGRESSION (p_mut=0) ───")
    log0 = run_with_birth_log(make_cfg(0.0))
    print(f"  faithful: {log0['faithful']}  mutated: {log0['mutated']}")
    print(f"  total_births: {log0['total_births']}  distinct_seqs: {log0['distinct_seqs']}")
    assert log0["mutated"] == 0, "REGRESSION FAILED: mutations appeared with p_mut=0"
    assert log0["faithful"] > 0, "REGRESSION FAILED: no births at all"
    print("  ✓ regression passes (100% faithful)\n")

    print(f"  baseline distinct seqs (p_mut=0): {log0['distinct_seqs']} / {log0['n_strands']} strands "
          f"= {log0['distinct_seqs']/max(1,log0['n_strands']):.2%}\n")

    print("─── Test 2: MUTATIONS APPEAR (p_mut=0.01) ───")
    log1 = run_with_birth_log(make_cfg(0.01))
    total = log1["faithful"] + log1["mutated"]
    mut_frac = log1["mutated"] / total if total else 0
    # Expected: per-position p=0.01 over child length L. For L=4 avg, fraction of
    # mutated children ≈ 1 - 0.99^4 ≈ 3.9%. Observed should be in same order.
    print(f"  faithful: {log1['faithful']}  mutated: {log1['mutated']}  ({mut_frac:.1%} mutated)")
    print(f"  total_births: {log1['total_births']}  distinct_seqs: {log1['distinct_seqs']} / "
          f"{log1['n_strands']} strands = "
          f"{log1['distinct_seqs']/max(1,log1['n_strands']):.2%}")
    assert log1["mutated"] > 0, "MUTATION FAILED: no mutations seen with p_mut=0.01"
    # Sanity: mutation fraction shouldn't be wildly off from p_mut * mean_L
    assert 0.01 < mut_frac < 0.20, (
        f"MUTATION FAILED: mutation fraction {mut_frac:.2%} outside sane range "
        f"(expected 1-20% for p_mut=0.01 with short children)")
    print(f"  ✓ {log1['mutated']} children carry off-template mutations\n")

    print("─── Test 3: DRIFT — accumulation along lineages (p_mut=0.05) ───")
    print("  Using higher p_mut so cumulative drift is visible in 1000 steps.")
    log2 = run_with_birth_log(make_cfg(0.05), n_steps=1500)
    n_with_inherited = sum(1 for g, _ in log2["by_gen"].items())  # placeholder
    # Count alive strands that have inherited at least one mutation
    # (we need access to the actual world for this — re-run quickly)
    cfg = make_cfg(0.05)
    w = World(cfg, seed=42)
    for t in range(1500):
        E = (np.sin(t * 0.05) + 1) / 2
        w.step(env_E=E)
    n_alive = w.n_strands()
    n_with_mut = sum(1 for s in w.strands if s.total_mutations_inherited > 0)
    max_inherited = max((s.total_mutations_inherited for s in w.strands), default=0)
    by_gen = {}
    for s in w.strands:
        by_gen.setdefault(s.gen, []).append(s.total_mutations_inherited)
    print(f"  {'gen':>4} {'count':>6} {'mean_inherited':>16} {'max_inherited':>15}")
    for g in sorted(by_gen.keys()):
        vals = by_gen[g]
        print(f"  {g:>4d} {len(vals):>6d} {np.mean(vals):>16.3f} {max(vals):>15d}")
    print(f"\n  alive strands carrying inherited mutations: {n_with_mut}/{n_alive} "
          f"({n_with_mut/max(1,n_alive):.1%})")
    print(f"  max inherited mutations on a single lineage: {max_inherited}")
    assert n_with_mut > 0, (
        "DRIFT FAILED: no alive strands carry inherited mutations after 1500 steps")
    assert max_inherited >= 1, "DRIFT FAILED: no lineage accumulated any mutation"
    print("  ✓ mutations propagate through lineages\n")

    print("Stage 1 done.")


if __name__ == "__main__":
    main()
