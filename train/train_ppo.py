"""PPO training on BIOTRON v2 (catalysis + lipids + membranes).

Run as a module so the env package import works:
    python3 -m train.train_ppo --smoke      # quick 5K-timestep sanity test
    python3 -m train.train_ppo              # full 500K-timestep training

The agent controls one continuous action: env_E ∈ [0, 1].
Reward is compartment-count based (see env/biotron_env.py for details).
Episode is fixed-length (default 2000 steps).

Outputs (in models/):
    biotron_ppo_v2_<TIMESTAMP>.zip   — saved policy
"""
from __future__ import annotations
import argparse
import os
import time
from pathlib import Path

import numpy as np
from stable_baselines3 import PPO
from stable_baselines3.common.vec_env import SubprocVecEnv, DummyVecEnv
from stable_baselines3.common.callbacks import BaseCallback
from stable_baselines3.common.monitor import Monitor

from env import BiotronEnv
from env.biotron_env import make_default_v2_cfg
from sim import WorldConfig


def make_env_factory(cfg: WorldConfig, episode_steps: int):
    def _thunk():
        return Monitor(BiotronEnv(world_cfg=cfg, episode_steps=episode_steps))
    return _thunk


class CompartmentLogger(BaseCallback):
    """Periodically print mean episode reward (≈ cumulative compartment count) and
    a rough estimate of the current policy's mean action (= preferred temperature)."""
    def __init__(self):
        super().__init__()

    def _on_step(self) -> bool:
        return True

    def _on_rollout_end(self) -> None:
        ep_info = self.model.ep_info_buffer
        if not ep_info:
            return
        rewards = [ep["r"] for ep in ep_info]
        lens = [ep["l"] for ep in ep_info]
        mean_r = float(np.mean(rewards))
        mean_l = float(np.mean(lens))
        n = self.num_timesteps
        try:
            act, _ = self.model.predict(self.locals["obs_tensor"].cpu().numpy(), deterministic=True)
            mean_action = float(np.mean(act))
            mean_E = (mean_action + 1.0) / 2.0  # rescale [-1,1] → [0,1]
        except Exception:
            mean_E = float("nan")
        print(f"  step={n:>9d}  mean_ep_reward={mean_r:>7.1f}  "
              f"ep_len={mean_l:>5.0f}  mean_E≈{mean_E:.2f}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--timesteps", type=int, default=500_000)
    parser.add_argument("--n-envs", type=int, default=8)
    parser.add_argument("--episode-steps", type=int, default=2000)
    parser.add_argument("--smoke", action="store_true",
                        help="Quick smoke test: 4000 timesteps, 2 envs, 500-step episodes")
    args = parser.parse_args()

    if args.smoke:
        args.timesteps = 4_000
        args.n_envs = 2
        args.episode_steps = 500

    cfg = make_default_v2_cfg()
    print(f"BIOTRON v2 RL — grid {cfg.grid_h}×{cfg.grid_w}, "
          f"k_lipid={cfg.k_lipid}, k_poly={cfg.k_poly_boost}, "
          f"chemo={cfg.lipid_chemotaxis}, seal={cfg.membrane_seal}", flush=True)
    print(f"Training: {args.timesteps:,} timesteps, "
          f"{args.n_envs} parallel envs, episode={args.episode_steps} steps", flush=True)

    factory = make_env_factory(cfg, args.episode_steps)
    if args.smoke or args.n_envs == 1:
        vec = DummyVecEnv([factory for _ in range(args.n_envs)])
    else:
        vec = SubprocVecEnv([factory for _ in range(args.n_envs)])

    model = PPO(
        "MlpPolicy",
        vec,
        verbose=0,
        n_steps=512,
        batch_size=128,
        n_epochs=10,
        gamma=0.997,         # very long-horizon — compartments take many steps to form
        gae_lambda=0.95,
        clip_range=0.2,
        learning_rate=3e-4,
        ent_coef=0.02,       # bump exploration since reward signal is sparse early
        policy_kwargs=dict(net_arch=[64, 64]),
        device="cpu",
    )

    out_dir = Path(__file__).parent.parent / "models"
    out_dir.mkdir(exist_ok=True)
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    save_name = "biotron_ppo_v2_smoke" if args.smoke else f"biotron_ppo_v2_{timestamp}"
    save_path = out_dir / save_name

    callback = CompartmentLogger()
    t0 = time.time()
    model.learn(total_timesteps=args.timesteps, callback=callback, progress_bar=False)
    dt = time.time() - t0
    print(f"\nTraining done in {dt:.0f}s ({args.timesteps/dt:.0f} steps/s)", flush=True)

    model.save(str(save_path))
    print(f"Saved → {save_path}.zip", flush=True)
    vec.close()


if __name__ == "__main__":
    main()
