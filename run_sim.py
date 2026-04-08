"""Smoke test for the BIOTRON Python simulation.

Runs N steps with a sin-wave energy schedule (the v5d analog) and prints
periodic stats. No RL, no UI — just verifies the physics doesn't crash and
produces qualitatively reasonable behavior.
"""
import numpy as np
from sim import World, WorldConfig

cfg = WorldConfig(
    grid_h=32,
    grid_w=32,
    initial_mono_per_cell=15,
    diffuse_p=0.10,
    spont_rate=0.005,
    hbond_base_rate=0.6,
    stacking_bonus=3.0,
    poly_base_rate=0.05,
    min_birth_length=3,
    max_strand_length=24,
    max_strands=4000,
)

world = World(cfg, seed=42)
print(f"Initial:  free={world.total_free():,}  strands={world.n_strands()}")
print(f"Grid: {cfg.grid_h}×{cfg.grid_w}, energy schedule: sin wave")
print("-" * 76)
print(f"{'step':>5} {'E':>5} {'free':>7} {'strands':>8} {'births':>7} {'avg_L':>6} {'max_L':>6}")

N_STEPS = 500
for step in range(1, N_STEPS + 1):
    E = (np.sin(step * 0.05) + 1) / 2
    world.step(env_E=E)

    if step % 25 == 0:
        n = world.n_strands()
        if n > 0:
            lens = [s.length for s in world.strands]
            avg_len = np.mean(lens)
            max_len = max(lens)
        else:
            avg_len = 0
            max_len = 0
        print(f"{step:>5d} {E:>5.2f} {world.total_free():>7,d} {n:>8d} "
              f"{world.total_births:>7d} {avg_len:>6.1f} {max_len:>6d}")

print("-" * 76)
print(f"Final births: {world.total_births}")
print(f"Final strands: {world.n_strands()}")
if world.strands:
    longest = max(world.strands, key=lambda s: s.length)
    print(f"Longest: {longest!r}")
    repls = [s for s in world.strands if s.repl_count > 0]
    if repls:
        print(f"Replicators: {len(repls)}")
        for s in sorted(repls, key=lambda s: -s.repl_count)[:5]:
            print(f"  ×{s.repl_count} {s!r}")
