"""Evaluate the learned PPO policy vs hand-tuned baselines on the same env.

Runs each policy for 5 episodes × 2000 steps and reports mean rewards +
additional metrics (big compartments, max observed, catalyst enrichment).
"""
import sys
import numpy as np
from stable_baselines3 import PPO
from env import BiotronEnv
from env.biotron_env import make_default_v2_cfg
from analysis.compartments import find_compartments
from analysis.compartment_lineages import strands_in_compartment


def run_one(policy_fn, n_episodes=5, episode_steps=2000, seed_base=100):
    rewards = []
    max_bigs = []
    peak_hammings = []
    for ep in range(n_episodes):
        env = BiotronEnv(episode_steps=episode_steps)
        obs, _ = env.reset(seed=seed_base + ep)
        total_r = 0.0
        max_big = 0
        done = False
        while not done:
            action = policy_fn(obs)
            obs, r, term, trunc, info = env.step(action)
            total_r += r
            if info["n_big_compartments"] > max_big:
                max_big = info["n_big_compartments"]
            done = term or trunc
        rewards.append(total_r)
        max_bigs.append(max_big)
    return {
        "mean_reward": float(np.mean(rewards)),
        "std_reward": float(np.std(rewards)),
        "max_big": max(max_bigs),
        "mean_max_big": float(np.mean(max_bigs)),
    }


def const_policy(E_val):
    a = 2 * E_val - 1  # rescale [0,1] → [-1,1]
    return lambda obs: np.array([a], dtype=np.float32)


def random_policy(obs):
    return np.array([np.random.uniform(-1, 1)], dtype=np.float32)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 eval_v2_policies.py <model_path.zip>")
        sys.exit(1)
    model_path = sys.argv[1]
    print(f"Loading: {model_path}")
    model = PPO.load(model_path)

    def agent_policy(obs):
        action, _ = model.predict(obs, deterministic=True)
        return action

    policies = {
        "agent":        agent_policy,
        "const_E000":   const_policy(0.00),
        "const_E025":   const_policy(0.25),
        "const_E050":   const_policy(0.50),
        "random":       random_policy,
    }

    print(f"\n{'policy':>14} {'mean_reward':>13} {'std':>8} {'max_big':>9} {'avg_max_big':>12}")
    print("-" * 65)
    for name, fn in policies.items():
        r = run_one(fn, n_episodes=8, episode_steps=2000)
        print(f"{name:>14} {r['mean_reward']:>13.0f} {r['std_reward']:>8.0f} "
              f"{r['max_big']:>9d} {r['mean_max_big']:>12.1f}", flush=True)


if __name__ == "__main__":
    main()
