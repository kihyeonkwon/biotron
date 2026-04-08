"""Time-series of population dynamics under fixed energy schedules.

Goal: confirm that the births we see in the sweep are part of a real
exponential growth phase (logistic-curve-like), not just random bursts.

Tracks per-step:
  - n_strands (population)
  - births_per_step (rate)
  - n_long (length ≥ 4 — "interesting" replicators)
  - max_gen (lineage depth)

Outputs ASCII plots of each curve and a doubling-time estimate for the
exponential phase.
"""
import numpy as np
from sim import World, WorldConfig

# ─── Reuse schedule helpers ───
def sched_constant(E_val):
    return lambda step: E_val
def sched_sin(period, lo=0.0, hi=1.0):
    return lambda step: lo + (hi - lo) * (np.sin(2 * np.pi * step / period) + 1) / 2

SCHEDULES = {
    "constant_03": sched_constant(0.30),
    "sin_mild":    sched_sin(period=200, lo=0.15, hi=0.75),
    "constant_05": sched_constant(0.50),
}


def run_timeseries(name, sched_fn, n_steps=2000, seed=42):
    cfg = WorldConfig(
        grid_h=48, grid_w=48,
        initial_mono_per_cell=15,
        spont_rate=0.002,           # less explosive seeding
        hbond_base_rate=0.6,
        stacking_bonus=3.0,
        poly_base_rate=0.05,
        extend_rate=0.02,
        min_birth_length=3,
        max_strand_length=24,
        max_strands=12000,          # higher cap to see saturation, not artificial truncation
    )
    w = World(cfg, seed=seed)
    n_strands = np.zeros(n_steps, dtype=np.int32)
    births    = np.zeros(n_steps, dtype=np.int32)
    n_long    = np.zeros(n_steps, dtype=np.int32)
    max_gen   = np.zeros(n_steps, dtype=np.int32)
    cum_births = 0
    for t in range(n_steps):
        E = sched_fn(t)
        w.step(env_E=E)
        n_strands[t] = w.n_strands()
        new_births = w.total_births - cum_births
        births[t] = new_births
        cum_births = w.total_births
        n_long[t] = sum(1 for s in w.strands if s.length >= 4)
        max_gen[t] = max((s.gen for s in w.strands), default=0)
    return {
        "name": name,
        "n_strands": n_strands,
        "births": births,
        "n_long": n_long,
        "max_gen": max_gen,
        "total_births": cum_births,
        "world": w,
    }


def ascii_plot(data, label, width=70, height=12, log=False):
    d = np.asarray(data, dtype=float)
    if log:
        d = np.log(np.maximum(d, 1))
    n = len(d)
    # Bin x-axis
    bin_size = max(1, n // width)
    binned = np.array([d[i:i+bin_size].mean() for i in range(0, n, bin_size)])[:width]
    lo, hi = float(binned.min()), float(binned.max())
    rng = hi - lo if hi > lo else 1.0
    grid = [[' '] * len(binned) for _ in range(height)]
    for i, v in enumerate(binned):
        h = int((v - lo) / rng * (height - 1))
        grid[height - 1 - h][i] = '█'
    print(f"\n  {label}  range=[{lo:.1f}, {hi:.1f}]" + ("  (log)" if log else ""))
    for row in grid:
        print("    " + "".join(row))
    print("    " + "─" * len(binned))
    print(f"    step 0{' ' * (len(binned) - 12)}step {n}")


def doubling_time(series, window=50):
    """Estimate doubling time during the exponential phase.
    Find the window with maximal positive log-slope."""
    s = np.asarray(series, dtype=float)
    s = np.maximum(s, 1)
    log_s = np.log(s)
    n = len(log_s)
    if n < window + 1:
        return None
    best_slope = 0.0
    best_start = 0
    for i in range(0, n - window):
        slope = (log_s[i + window] - log_s[i]) / window
        if slope > best_slope:
            best_slope = slope
            best_start = i
    if best_slope <= 0:
        return None
    dt = np.log(2) / best_slope
    return dt, best_start, best_start + window


def main():
    import time
    print("Running 3 schedules × 2000 steps × 48×48 grid...\n", flush=True)
    results = {}
    for name, fn in SCHEDULES.items():
        t0 = time.time()
        results[name] = run_timeseries(name, fn)
        print(f"  {name}: total_births={results[name]['total_births']}, "
              f"final_strands={results[name]['n_strands'][-1]}, "
              f"final_max_gen={results[name]['max_gen'][-1]}, "
              f"({time.time()-t0:.1f}s)", flush=True)

    print("\n" + "=" * 78)
    for name, r in results.items():
        print(f"\n──── {name} ────")
        ascii_plot(r["n_strands"], "strand population (linear)", height=10)
        ascii_plot(r["n_strands"], "strand population (log)", height=10, log=True)
        # Smoothed birth rate
        births_smooth = np.convolve(r["births"], np.ones(20)/20, mode="same")
        ascii_plot(births_smooth, "births/step (20-step moving avg)", height=8)
        ascii_plot(r["max_gen"], "max generation alive", height=6)
        dt = doubling_time(r["n_strands"], window=50)
        if dt:
            d, a, b = dt
            print(f"\n  Best 50-step exponential window: steps {a}–{b},  doubling time ≈ {d:.0f} steps")
        else:
            print("\n  No exponential growth window found.")
        # Sustained metric: how many consecutive steps had births > 0
        nonzero = r["births"] > 0
        max_run = 0
        run = 0
        for v in nonzero:
            if v:
                run += 1
                max_run = max(max_run, run)
            else:
                run = 0
        print(f"  Longest consecutive-birth streak: {max_run} steps")


if __name__ == "__main__":
    main()
