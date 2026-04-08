"""Evaluate fixed-policy baselines using the same env as RL training.
This gives us numbers to beat: any RL policy must outperform constant baselines."""
import numpy as np
from env import BiotronEnv

EPISODE_STEPS = 1500
N_EPISODES = 5

POLICIES = {
    "const_E020": lambda obs: np.array([2 * 0.20 - 1], dtype=np.float32),
    "const_E030": lambda obs: np.array([2 * 0.30 - 1], dtype=np.float32),
    "const_E050": lambda obs: np.array([2 * 0.50 - 1], dtype=np.float32),
    "const_E070": lambda obs: np.array([2 * 0.70 - 1], dtype=np.float32),
    "pcr_60c_12h": lambda obs, _state={"t": 0}: (
        np.array([2 * (0.95 if (_state["t"] % 72) >= 60 else 0.20) - 1], dtype=np.float32),
        _state.update({"t": _state["t"] + 1}),
    )[0],
    "random": lambda obs: np.array([np.random.uniform(-1, 1)], dtype=np.float32),
}


def run_policy(name, policy_fn, n_eps=N_EPISODES, episode_steps=EPISODE_STEPS):
    env = BiotronEnv(episode_steps=episode_steps)
    rewards = []
    for ep in range(n_eps):
        # reset internal state for stateful policies
        if hasattr(policy_fn, "__defaults__"):
            for d in (policy_fn.__defaults__ or ()):
                if isinstance(d, dict) and "t" in d:
                    d["t"] = 0
        obs, _ = env.reset(seed=42 + ep)
        total = 0.0
        done = False
        while not done:
            action = policy_fn(obs)
            obs, r, term, trunc, _ = env.step(action)
            total += r
            done = term or trunc
        rewards.append(total)
    return float(np.mean(rewards)), float(np.std(rewards))


def main():
    print(f"Each policy: {N_EPISODES} episodes × {EPISODE_STEPS} steps")
    print(f"{'policy':>14}  {'mean_births':>12}  {'± std':>8}")
    print("-" * 42)
    for name, fn in POLICIES.items():
        m, s = run_policy(name, fn)
        print(f"{name:>14}  {m:>12.0f}  {s:>8.0f}", flush=True)


if __name__ == "__main__":
    main()
