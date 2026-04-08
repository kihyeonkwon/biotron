"""Compartment lineage analysis (Stage 6).

Computes per-compartment statistics:
  - resident strand population
  - consensus sequence (most common monomer at each position, length-3 windows)
  - mean motif counts
  - pairwise Hamming distance between compartments

Used by the Stage 6 acceptance test.
"""
from __future__ import annotations
import numpy as np
from collections import Counter
from sim import catalysis as cat_mod


def strands_in_compartment(world, interior: set) -> list:
    """Return strands whose (x,y) is inside the given interior cell set."""
    return [s for s in world.strands if (s.x, s.y) in interior]


def consensus_sequence(strands: list, max_len: int = 8) -> tuple:
    """Pad/truncate strands to a fixed length and return per-position majority."""
    if not strands:
        return tuple()
    L = max_len
    counts = [Counter() for _ in range(L)]
    for s in strands:
        m = s.mono.tolist()
        for i in range(min(L, len(m))):
            counts[i][m[i]] += 1
    return tuple(c.most_common(1)[0][0] if c else -1 for c in counts)


def hamming(a: tuple, b: tuple) -> int:
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] != b[i] and a[i] != -1 and b[i] != -1)


def compartment_stats(world, interiors: list) -> list:
    """Per-compartment summary list."""
    stats = []
    for i, comp in enumerate(interiors):
        strands = strands_in_compartment(world, comp)
        if not strands:
            continue
        cons = consensus_sequence(strands)
        mean_lipid_motif = float(np.mean([s._cat_lipid for s in strands]))
        mean_poly_motif  = float(np.mean([s._cat_poly for s in strands]))
        stats.append({
            "id": i,
            "n_cells": len(comp),
            "n_strands": len(strands),
            "consensus": cons,
            "mean_cat_lipid": mean_lipid_motif,
            "mean_cat_poly": mean_poly_motif,
        })
    return stats


def pairwise_hamming(stats: list) -> list:
    """Return list of (i, j, hamming_distance) tuples."""
    out = []
    for i in range(len(stats)):
        for j in range(i + 1, len(stats)):
            d = hamming(stats[i]["consensus"], stats[j]["consensus"])
            out.append((stats[i]["id"], stats[j]["id"], d))
    return out
