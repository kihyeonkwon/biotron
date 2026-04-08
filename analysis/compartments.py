"""Compartment finder: flood-fill on non-membrane cells.

A compartment = a connected component of non-membrane cells that does NOT
touch the grid boundary. The component touching the grid boundary is the
"exterior" (the rest of the soup).

Used by Stage 4-6 verification to detect emergent proto-cell interiors.
"""
from __future__ import annotations
import numpy as np


def find_compartments(is_membrane: np.ndarray) -> tuple[set[tuple[int, int]], list[set[tuple[int, int]]]]:
    """Run flood-fill on the non-membrane region and return (exterior, interiors).

    exterior: set of (x,y) cells in the connected component touching the grid edge.
              May be empty if the entire boundary is membrane (degenerate case).
    interiors: list of cell-sets, one per enclosed (non-boundary-touching) component.
               Empty if no compartments exist.
    """
    H, W = is_membrane.shape
    visited = np.zeros_like(is_membrane, dtype=bool)
    components: list[set[tuple[int, int]]] = []
    boundary_components: list[int] = []  # indices into components

    for i in range(H):
        for j in range(W):
            if is_membrane[i, j] or visited[i, j]:
                continue
            # BFS
            comp: set[tuple[int, int]] = set()
            stack = [(i, j)]
            touches_boundary = False
            while stack:
                x, y = stack.pop()
                if visited[x, y] or is_membrane[x, y]:
                    continue
                visited[x, y] = True
                comp.add((x, y))
                if x == 0 or x == H - 1 or y == 0 or y == W - 1:
                    touches_boundary = True
                # 4-neighbors
                if x > 0 and not visited[x - 1, y] and not is_membrane[x - 1, y]:
                    stack.append((x - 1, y))
                if x < H - 1 and not visited[x + 1, y] and not is_membrane[x + 1, y]:
                    stack.append((x + 1, y))
                if y > 0 and not visited[x, y - 1] and not is_membrane[x, y - 1]:
                    stack.append((x, y - 1))
                if y < W - 1 and not visited[x, y + 1] and not is_membrane[x, y + 1]:
                    stack.append((x, y + 1))
            components.append(comp)
            if touches_boundary:
                boundary_components.append(len(components) - 1)

    # Exterior = union of all boundary-touching components
    exterior: set[tuple[int, int]] = set()
    for idx in boundary_components:
        exterior |= components[idx]
    interiors = [c for i, c in enumerate(components) if i not in boundary_components]
    return exterior, interiors


def compartment_stats(exterior: set, interiors: list[set]) -> dict:
    """Summary statistics for diagnostics."""
    return {
        "n_compartments": len(interiors),
        "exterior_cells": len(exterior),
        "interior_total_cells": sum(len(c) for c in interiors),
        "largest_interior": max((len(c) for c in interiors), default=0),
        "interior_sizes": sorted((len(c) for c in interiors), reverse=True),
    }
