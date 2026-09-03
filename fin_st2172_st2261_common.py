#!/usr/bin/env python3
"""Shared finite mathematics for FIN ST2172--ST2261.

The helpers distinguish a reducible three-cycle statistic of the strict pair
kernel from genuinely irreducible third-order state information.  No physical
interpretation, selector, dimensional scale, or legacy role is assumed.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any
import numpy as np


ROOT = Path(__file__).resolve().parent
N = 12
VERTICES = tuple(range(N))
EDGES = tuple(itertools.combinations(VERTICES, 2))
TRIANGLES = tuple(itertools.combinations(VERTICES, 3))
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}


def strict_A_W() -> tuple[np.ndarray, np.ndarray]:
    """Return the strict Dirichlet generator and its positive off-diagonal W."""
    omega, phi, eta = 0.18575, 0.16250, 9 / 5
    weights = {d: math.cos(omega * d + phi) / (1 + d**eta) for d in range(1, 7)}
    s = 2 * sum(weights[d] for d in range(1, 6)) + weights[6]
    A = np.array([[s if i == j else -weights[min((i-j) % N, (j-i) % N)]
                   for j in range(N)] for i in range(N)])
    W = -A.copy()
    np.fill_diagonal(W, 0.0)
    return A, W


def write_packet(program: int, name: str, status: str, boundary: str,
                 payload: dict[str, Any]) -> dict[str, Any]:
    path = ROOT / f"FIN_ST{program}_{name}.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return {"program": f"ST{program}", "object": name,
            "packet_file": path.name,
            "packet_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            **payload, "status": status, "boundary": boundary}


def write_round(lo: int, hi: int, results: dict[str, Any]) -> None:
    assert hi - lo + 1 == 15
    assert list(results) == [f"ST{k}" for k in range(lo, hi + 1)]
    (ROOT / f"FIN_ST{lo}_ST{hi}_Results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n")
    with (ROOT / f"FIN_ST{lo}_ST{hi}_Summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["program", "status", "object", "boundary"])
        for key, value in results.items():
            writer.writerow([key, value["status"], value["object"], value["boundary"]])


def d12_permutations() -> list[tuple[str, int, tuple[int, ...]]]:
    out = []
    for k in VERTICES:
        out.append(("rotation", k, tuple((i + k) % N for i in VERTICES)))
        out.append(("reflection", k, tuple((k - i) % N for i in VERTICES)))
    return out


def permutation_matrix(p: tuple[int, ...]) -> np.ndarray:
    P = np.zeros((N, N))
    for old, new in enumerate(p):
        P[new, old] = 1.0
    return P


def cochains() -> tuple[np.ndarray, np.ndarray]:
    d0 = np.zeros((len(EDGES), N))
    d1 = np.zeros((len(TRIANGLES), len(EDGES)))
    for r, (i, j) in enumerate(EDGES):
        d0[r, i], d0[r, j] = -1.0, 1.0
    for r, (i, j, k) in enumerate(TRIANGLES):
        d1[r, EDGE_INDEX[(i, j)]] = 1.0
        d1[r, EDGE_INDEX[(j, k)]] = 1.0
        d1[r, EDGE_INDEX[(i, k)]] = -1.0
    return d0, d1


def triangle_orbits() -> list[set[tuple[int, int, int]]]:
    remaining = set(TRIANGLES)
    orbits = []
    while remaining:
        tri = min(remaining)
        orbit = set()
        for _, _, p in d12_permutations():
            orbit.add(tuple(sorted(p[i] for i in tri)))
        orbit &= set(TRIANGLES)
        orbits.append(orbit)
        remaining -= orbit
    return orbits


def cycle_field(W: np.ndarray) -> np.ndarray:
    return np.array([W[i, j] * W[j, k] * W[k, i] for i, j, k in TRIANGLES])


def normalized_positive_power(tau: np.ndarray, exponent: float) -> np.ndarray:
    h = np.power(tau, exponent)
    return h / np.mean(h)


def hodge_spectrum(h: np.ndarray) -> np.ndarray:
    d0, d1 = cochains()
    L1 = d0 @ d0.T + d1.T @ np.diag(h) @ d1
    return np.linalg.eigvalsh(L1)


def fourier_state(k: int) -> np.ndarray:
    x = np.arange(N)
    return np.exp(2j * np.pi * k * x / N) / np.sqrt(N)


def edge_current(rho: np.ndarray, W: np.ndarray) -> float:
    return float(sum(2.0 * np.imag(rho[i, (i + 1) % N] * W[(i + 1) % N, i]) for i in VERTICES))


def hidden_synergy_distribution(epsilon: float) -> dict[tuple[int, int, int], float]:
    if abs(epsilon) > 1.0:
        raise ValueError("positivity requires |epsilon| <= 1")
    return {
        x: (1.0 + epsilon * x[0] * x[1] * x[2]) / 8.0
        for x in itertools.product((-1, 1), repeat=3)
    }


def moments(distribution: dict[tuple[int, int, int], float]) -> dict[str, object]:
    means = [sum(p * x[i] for x, p in distribution.items()) for i in range(3)]
    pairs = [sum(p * x[i] * x[j] for x, p in distribution.items()) for i, j in ((0, 1), (0, 2), (1, 2))]
    triple = sum(p * x[0] * x[1] * x[2] for x, p in distribution.items())
    return {"means": means, "pairs": pairs, "triple": triple, "mass": sum(distribution.values()),
            "minimum_probability": min(distribution.values())}
