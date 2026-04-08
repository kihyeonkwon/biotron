"""Core physics constants and Boltzmann mechanics for BIOTRON.

5 physics rules (these are the ONLY rules in the simulation):
  1. H-bond  (Donor ↔ Acceptor)   — weak, reversible
  2. Covalent bond                — strong, *weakly* reversible
  3. Stacking                     — adjacent bonded positions stabilize each other
  4. Energy                       — environment thermal energy (RL agent controls)
  5. Boltzmann survival           — every break is P = exp(-bondE / kT)

Template-directed replication is NOT coded. It must emerge.
"""
import numpy as np

# ── Bond energies (arbitrary units; what matters is their ratios to kT) ──
E_COV_PARENT   = 1.20   # parent backbone covalent (mature, more stable)
E_COV_DAUGHTER = 0.95   # daughter backbone covalent (newly polymerized — slightly weaker)
E_HB           = 0.30   # H-bond
E_STACK        = 0.15   # stacking bonus per adjacent bonded neighbor

# ── Thermal energy as a function of environment energy ∈ [0, 1] ──
KT_BASE  = 0.05
KT_RANGE = 0.30

def kT(env_E):
    """Thermal energy. env_E in [0, 1]."""
    return KT_BASE + KT_RANGE * env_E

def break_prob(bond_E, env_E):
    """Boltzmann break probability per step. P(break) = exp(-bondE / kT)."""
    return float(np.exp(-bond_E / kT(env_E)))

# ── Monomer types ──
# A monomer has 2 sites, each Donor(1) or Acceptor(0).
# 4 types total: encode as integer 0..3 = (s0 * 2 + s1)
#   0 = (0,0)   1 = (0,1)   2 = (1,0)   3 = (1,1)
# Complement (D↔A everywhere): bit-flip both → 0↔3, 1↔2
N_TYPES = 4
COMPLEMENT = np.array([3, 2, 1, 0], dtype=np.int8)

# Human-readable labels for debugging only (no semantic meaning to the simulation)
LABELS = {0: "α", 1: "β", 2: "β\u0304", 3: "α\u0304"}
