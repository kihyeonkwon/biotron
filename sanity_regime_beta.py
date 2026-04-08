"""Regime β sanity check: prove cycling is REQUIRED when H-bonds are strong.

Hypothesis:
  - Regime α (E_HB=0.30, current): constant warm works → cycling not needed.
  - Regime β (E_HB=0.80): constant warm fails (daughters can't lift off) →
    only schedules with hot pulses succeed.

If this contrast shows up cleanly, the simulation is responsive to schedule
choice and RL on Regime β should later be able to discover PCR-like cycling.
"""
import numpy as np
from sim import World, WorldConfig

# ─── Regimes ───
def make_regime_alpha(**overrides):
    """Default physics — weak H-bonds. Constant warm works."""
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
        # default energies (Regime α)
        e_hb=0.30,
        e_stack=0.15,
        e_cov_parent=1.20,
        e_cov_daughter=0.95,
        **overrides,
    )

def make_regime_beta(**overrides):
    """Strong H-bond physics — daughters stuck at warm, hot pulse required."""
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
        e_hb=0.80,           # ~2.7× stronger H-bond
        e_stack=0.20,        # slightly stronger stacking too
        e_cov_parent=1.80,   # parent more durable, can survive hot pulses
        e_cov_daughter=1.50, # daughter cov also more durable
    )
    base.update(overrides)
    return WorldConfig(**base)


# ─── Schedules ───
def sched_const(E):
    return lambda step: E

def sched_pcr_square(cold_E, hot_E, cold_steps, hot_steps):
    cycle = cold_steps + hot_steps
    return lambda step: hot_E if (step % cycle) >= cold_steps else cold_E

SCHEDULES = {
    "const_warm_E020":    sched_const(0.20),
    "const_warm_E030":    sched_const(0.30),
    "const_warm_E050":    sched_const(0.50),
    "pcr_100c_8h":        sched_pcr_square(0.20, 0.95, 100, 8),
    "pcr_60c_12h":        sched_pcr_square(0.20, 0.95, 60, 12),
    "pcr_40c_6h":         sched_pcr_square(0.20, 0.95, 40, 6),
}


def run_one(regime_name, regime_fn, sched_name, sched_fn, n_steps=1500, seed=42):
    cfg = regime_fn()
    w = World(cfg, seed=seed)
    for t in range(n_steps):
        w.step(env_E=sched_fn(t))
    return {
        "regime": regime_name,
        "schedule": sched_name,
        "births": w.total_births,
        "strands": w.n_strands(),
        "max_gen": max((s.gen for s in w.strands), default=0),
        "n_gen2": sum(1 for s in w.strands if s.gen >= 2),
        "n_long": sum(1 for s in w.strands if s.length >= 4),
    }


def main():
    import time
    print(f"\n{'regime':>10} {'schedule':>18} {'births':>8} {'strands':>8} "
          f"{'maxGen':>7} {'gen≥2':>7} {'len≥4':>7}", flush=True)
    print("-" * 70, flush=True)
    regimes = [("alpha", make_regime_alpha), ("beta", make_regime_beta)]
    for r_name, r_fn in regimes:
        for s_name, s_fn in SCHEDULES.items():
            t0 = time.time()
            r = run_one(r_name, r_fn, s_name, s_fn, n_steps=1500)
            print(f"{r['regime']:>10} {r['schedule']:>18} {r['births']:>8d} "
                  f"{r['strands']:>8d} {r['max_gen']:>7d} {r['n_gen2']:>7d} "
                  f"{r['n_long']:>7d}  ({time.time()-t0:.0f}s)", flush=True)
        print("-" * 70, flush=True)


if __name__ == "__main__":
    main()
