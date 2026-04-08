"""Sweep energy schedules to find one that produces sustained chain replication.

For each schedule, run a 1500-step simulation and report:
  - total_births
  - max_repl_count   (best single replicator)
  - n_replicators    (strands with repl_count > 0)
  - max_gen          (deepest generation reached — gen >= 2 means a child replicated)
  - n_gen2_plus      (strands with gen >= 2 — proves lineage depth)
  - n_gen3_plus      (gen >= 3 — proves chain of 3 replications)
"""
import numpy as np
from sim import World, WorldConfig

# ─────────── Schedules ───────────

def sched_constant(E_val):
    return lambda step: E_val

def sched_sin(period, lo=0.0, hi=1.0):
    return lambda step: lo + (hi - lo) * (np.sin(2 * np.pi * step / period) + 1) / 2

def sched_square(cold_E, hot_E, cold_steps, hot_steps):
    cycle = cold_steps + hot_steps
    def f(step):
        return hot_E if (step % cycle) >= cold_steps else cold_E
    return f

def sched_ramp(period, lo=0.05, hi=0.95):
    half = period / 2
    def f(step):
        s = step % period
        if s < half:
            return lo + (hi - lo) * (s / half)
        else:
            return hi - (hi - lo) * ((s - half) / half)
    return f


SCHEDULES = {
    "constant_03":      sched_constant(0.30),
    "constant_05":      sched_constant(0.50),
    "sin_fast":         sched_sin(period=126, lo=0.0, hi=1.0),
    "sin_slow":         sched_sin(period=300, lo=0.0, hi=1.0),
    "sin_mild":         sched_sin(period=200, lo=0.15, hi=0.75),
    "square_50_5":      sched_square(cold_E=0.15, hot_E=0.85, cold_steps=50, hot_steps=5),
    "square_80_5":      sched_square(cold_E=0.15, hot_E=0.85, cold_steps=80, hot_steps=5),
    "square_100_3":     sched_square(cold_E=0.15, hot_E=0.85, cold_steps=100, hot_steps=3),
    "square_30_8":      sched_square(cold_E=0.20, hot_E=0.80, cold_steps=30, hot_steps=8),
    "ramp_200":         sched_ramp(period=200, lo=0.05, hi=0.90),
}


def run_one(name, schedule_fn, n_steps=1000, seed=42):
    cfg = WorldConfig(
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
    )
    w = World(cfg, seed=seed)
    for step in range(n_steps):
        E = schedule_fn(step)
        w.step(env_E=E)

    strands = w.strands
    repls = [s for s in strands if s.repl_count > 0]
    max_repl = max((s.repl_count for s in strands), default=0)
    max_gen = max((s.gen for s in strands), default=0)
    n_gen2 = sum(1 for s in strands if s.gen >= 2)
    n_gen3 = sum(1 for s in strands if s.gen >= 3)

    # Lineage depth = max generation observed AT BIRTH (not just current)
    # We need history. Approximate by current strands' max gen.
    # For a more accurate measure, track births separately.
    return {
        "name": name,
        "total_births": w.total_births,
        "n_strands": len(strands),
        "n_replicators": len(repls),
        "max_repl": max_repl,
        "max_gen_alive": max_gen,
        "n_gen2_alive": n_gen2,
        "n_gen3_alive": n_gen3,
    }


def main():
    import time
    print(f"{'schedule':>16} {'births':>8} {'strands':>8} {'repls':>7} "
          f"{'maxRepl':>8} {'maxGen':>7} {'gen≥2':>7} {'gen≥3':>7} {'time':>6}", flush=True)
    print("-" * 84, flush=True)
    results = []
    for name, fn in SCHEDULES.items():
        t0 = time.time()
        r = run_one(name, fn, n_steps=1000, seed=42)
        dt = time.time() - t0
        results.append(r)
        print(f"{r['name']:>16} {r['total_births']:>8d} {r['n_strands']:>8d} "
              f"{r['n_replicators']:>7d} {r['max_repl']:>8d} {r['max_gen_alive']:>7d} "
              f"{r['n_gen2_alive']:>7d} {r['n_gen3_alive']:>7d} {dt:>5.1f}s", flush=True)
    print("-" * 84, flush=True)
    # Best by max_gen, then by max_repl, then by total_births
    best = max(results, key=lambda r: (r["max_gen_alive"], r["max_repl"], r["total_births"]))
    print(f"\nBest by lineage depth: {best['name']}")
    print(f"  max_gen={best['max_gen_alive']}, max_repl={best['max_repl']}, "
          f"total_births={best['total_births']}")


if __name__ == "__main__":
    main()
