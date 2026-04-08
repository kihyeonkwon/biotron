"""Gymnasium environment wrapping the BIOTRON v2 world.

Action: continuous E ∈ [0, 1] (rescaled from [-1, 1] for Gaussian policy friendliness).
Reward: number of "big" compartments (≥5 interior cells) at this step.
        Small baseline reward proportional to membrane formation so the agent
        gets non-zero signal in the lag phase before compartments emerge.
Observation: 9-dim feature vector covering strand population, lipid pool,
             membrane state, compartment count, mean catalyst strengths, and
             current energy.
Episode: fixed length (default 2000 steps).
"""
from collections import deque
import numpy as np
import gymnasium as gym
from gymnasium import spaces

from sim import World, WorldConfig
from analysis.compartments import find_compartments


# ─── Default v2.1 cfg (the one that hit 4/4 in Stage 6) ───
def make_default_v2_cfg() -> WorldConfig:
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
        # Stage 1
        p_mut=0.005,
        # Stage 2
        k_poly_boost=2.0,
        # Stage 3
        k_lipid=1.0,
        lipid_diffuse_p=0.02,
        lipid_chemotaxis=0.10,
        e_lipid=1.0,
        lipid_decay_attempt_rate=1.0,
        # Stage 4
        membrane_threshold_high=8,
        membrane_threshold_low=4,
        # Stage 5
        membrane_seal=0.95,
    )


# Observation normalization scales (rough)
N_STRANDS_SCALE      = 4000.0
TOTAL_FREE_SCALE     = 30000.0
TOTAL_LIPID_SCALE    = 2000.0
N_MEMBRANE_SCALE     = 200.0
N_COMPARTMENTS_SCALE = 10.0
CAT_LIPID_SCALE      = 1.0
CAT_POLY_SCALE       = 1.0
RECENT_BIRTH_SCALE   = 50.0


# Reward shaping constants
COMPARTMENT_MIN_SIZE       = 5      # only count "big" compartments
BIG_COMPARTMENT_REWARD     = 1.0    # per big compartment per step
SMALL_COMPARTMENT_REWARD   = 0.05   # per small compartment per step
MEMBRANE_FORMATION_BONUS   = 0.001  # per membrane cell per step (early signal)


class BiotronEnv(gym.Env):
    """v2.1 BIOTRON environment.

    Reset → fresh World with default v2 cfg
    Step  → advance one physics step under env_E = (action+1)/2
    Reward: weighted count of compartments + tiny membrane bonus
    Truncated when step_count >= episode_steps
    """

    metadata = {"render_modes": []}

    def __init__(
        self,
        world_cfg: WorldConfig | None = None,
        episode_steps: int = 2000,
        recent_window: int = 20,
        compartment_eval_every: int = 1,  # how often to recompute compartments (1 = every step)
    ):
        super().__init__()
        self.cfg = world_cfg or make_default_v2_cfg()
        self.episode_steps = episode_steps
        self.recent_window = recent_window
        self.compartment_eval_every = compartment_eval_every

        # Action: env_E rescaled from [-1, 1]
        self.action_space = spaces.Box(low=-1.0, high=1.0, shape=(1,), dtype=np.float32)
        # 9 features (see _obs)
        self.observation_space = spaces.Box(
            low=np.zeros(9, dtype=np.float32),
            high=np.array([5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 1.0], dtype=np.float32),
            dtype=np.float32,
        )

        self.world: World | None = None
        self.steps = 0
        self.recent_births: deque = deque(maxlen=recent_window)
        # Cached compartment count from the last evaluation
        self._n_big_compartments = 0
        self._n_small_compartments = 0

    # ─────────── Gym interface ───────────

    def reset(self, seed=None, options=None):
        super().reset(seed=seed)
        self.world = World(self.cfg, seed=seed)
        self.steps = 0
        self.recent_births = deque([0] * self.recent_window, maxlen=self.recent_window)
        self._n_big_compartments = 0
        self._n_small_compartments = 0
        return self._obs(), {}

    def step(self, action):
        a = float(np.clip(action[0], -1.0, 1.0))
        E = (a + 1.0) / 2.0
        births_before = self.world.total_births
        self.world.step(env_E=E)
        self.steps += 1
        births_this_step = self.world.total_births - births_before
        self.recent_births.append(births_this_step)

        # Compartment evaluation (cache between calls if compartment_eval_every > 1)
        if self.steps % self.compartment_eval_every == 0:
            _, interiors = find_compartments(self.world.is_membrane)
            self._n_big_compartments = sum(1 for c in interiors if len(c) >= COMPARTMENT_MIN_SIZE)
            self._n_small_compartments = len(interiors) - self._n_big_compartments

        n_membrane = int(self.world.is_membrane.sum())
        reward = (
            self._n_big_compartments * BIG_COMPARTMENT_REWARD
            + self._n_small_compartments * SMALL_COMPARTMENT_REWARD
            + n_membrane * MEMBRANE_FORMATION_BONUS
        )

        terminated = False
        truncated = self.steps >= self.episode_steps
        info = {
            "n_strands": self.world.n_strands(),
            "n_big_compartments": self._n_big_compartments,
            "n_small_compartments": self._n_small_compartments,
            "n_membrane_cells": n_membrane,
            "total_lipid": int(self.world.lipid.sum()),
            "births_this_step": births_this_step,
        }
        return self._obs(), reward, terminated, truncated, info

    # ─────────── Observation ───────────

    def _obs(self):
        w = self.world
        n = w.n_strands()
        if n > 0:
            mean_cat_lipid = float(np.mean([s._cat_lipid for s in w.strands]))
            mean_cat_poly  = float(np.mean([s._cat_poly  for s in w.strands]))
        else:
            mean_cat_lipid = 0.0
            mean_cat_poly = 0.0
        recent = sum(self.recent_births) / max(1, len(self.recent_births))
        n_membrane = int(w.is_membrane.sum())
        n_compartments = self._n_big_compartments + self._n_small_compartments
        return np.array([
            n / N_STRANDS_SCALE,
            w.total_free() / TOTAL_FREE_SCALE,
            int(w.lipid.sum()) / TOTAL_LIPID_SCALE,
            n_membrane / N_MEMBRANE_SCALE,
            n_compartments / N_COMPARTMENTS_SCALE,
            mean_cat_lipid / CAT_LIPID_SCALE,
            mean_cat_poly / CAT_POLY_SCALE,
            recent / RECENT_BIRTH_SCALE,
            float(w.env_E),
        ], dtype=np.float32)
