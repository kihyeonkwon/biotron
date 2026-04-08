"""Sequence → catalytic function (Stage 2).

We define a small set of trigram "motifs". A strand's catalytic strength for
a given reaction is the number of (overlapping) occurrences of the corresponding
motif in its sequence.

This is the simplest sequence-based mapping that:
  - is heritable (motifs are part of the sequence — replication carries them)
  - is mutable (a single substitution can create or destroy a motif site)
  - has a smooth-enough fitness landscape (small steps in sequence ↔ small steps
    in catalytic strength)
  - is O(L) per strand per reaction (cheap)

Motifs are short (length 3) so a length-10 strand has ~8 trigram windows. With
4 monomer types there are 64 possible trigrams, so a random strand has ~0.125
hits per motif on average — non-trivial signal but plenty of room for selection.
"""
from __future__ import annotations
import numpy as np

# ─── Reaction registry ───
# Each entry: name → np.int8 trigram. Hand-picked so they don't all share
# common subsequences and have somewhat different "shapes" (palindromic,
# alternating, mixed).

# Reaction names (string keys). Used by World/Strand for lookup.
REACTIONS = ("lipid_synth", "polymerize", "extend")

# Default motifs (Stage 2 uses just polymerize; lipid_synth is wired in Stage 3).
DEFAULT_MOTIFS = {
    "lipid_synth": np.array([0, 3, 0], dtype=np.int8),  # palindromic α–ᾱ–α
    "polymerize":  np.array([1, 2, 1], dtype=np.int8),  # palindromic β–β̄–β
    "extend":      np.array([0, 1, 2], dtype=np.int8),  # mixed run
}


def catalytic_strength(mono: np.ndarray, motif: np.ndarray) -> int:
    """Count overlapping occurrences of `motif` in `mono`. O(L * len(motif))."""
    L = len(mono)
    M = len(motif)
    if L < M:
        return 0
    count = 0
    # Vectorized: build all sliding windows then compare row-wise
    # For trigrams this is faster than a Python loop with non-trivial L
    if M == 3:
        # Specialize the common case
        eq0 = (mono[:L - 2] == motif[0])
        eq1 = (mono[1:L - 1] == motif[1])
        eq2 = (mono[2:L] == motif[2])
        return int(np.count_nonzero(eq0 & eq1 & eq2))
    # Fallback for other motif lengths
    for i in range(L - M + 1):
        if (mono[i:i + M] == motif).all():
            count += 1
    return count


def all_strengths(mono: np.ndarray, motifs: dict | None = None) -> dict:
    """Compute strength for every reaction. Returns {reaction: int}."""
    motifs = motifs or DEFAULT_MOTIFS
    return {name: catalytic_strength(mono, m) for name, m in motifs.items()}
