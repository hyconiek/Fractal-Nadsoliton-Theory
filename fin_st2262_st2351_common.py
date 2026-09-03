#!/usr/bin/env python3
"""Shared calculations for FIN ST2262--ST2351."""
from __future__ import annotations

import math
import numpy as np

from fin_st2172_st2261_common import (
    EDGES, TRIANGLES, ROOT, cycle_field, strict_A_W, triangle_orbits,
    write_packet, write_round,
)


def orbit_cycle_vector() -> np.ndarray:
    _, W = strict_A_W()
    tau = cycle_field(W)
    return np.array([np.mean([tau[TRIANGLES.index(t)] for t in orbit])
                     for orbit in triangle_orbits()])


def mean_normalize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return x / np.mean(x)


def power_map(x: np.ndarray, p: float) -> np.ndarray:
    return mean_normalize(np.asarray(x, dtype=float) ** p)


def entropy_of_mean_one_weights(x: np.ndarray) -> float:
    probability = mean_normalize(x) / len(x)
    return float(-sum(v * math.log(v) for v in probability if v > 0))


def strict_distance_weights() -> np.ndarray:
    _, W = strict_A_W()
    return np.array([W[0, d] for d in range(1, 7)])


__all__ = [
    "EDGES", "TRIANGLES", "ROOT", "cycle_field", "strict_A_W",
    "triangle_orbits", "write_packet", "write_round", "orbit_cycle_vector",
    "mean_normalize", "power_map", "entropy_of_mean_one_weights",
    "strict_distance_weights",
]
