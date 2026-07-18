#!/usr/bin/env python3
"""Deterministic checks for the FIN wave/diffusion observer analysis.

Uses only NumPy.  All times are dimensionless.  The script does not import
repository verdicts; the QW2118 radial values are declared explicitly below.
"""

from __future__ import annotations

import json
import math
from typing import Callable

import numpy as np


N = 12
RADIAL = np.array(
    [
        0.0,
        0.469985673,
        0.192043552,
        0.091428614,
        0.047029169,
        0.024131223,
        0.011070817,
    ],
    dtype=float,
)


def cyclic_matrix(radial: np.ndarray) -> np.ndarray:
    out = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            d = abs(i - j)
            d = min(d, N - d)
            out[i, j] = radial[d]
    return out


W = cyclic_matrix(RADIAL)
ROW_SUM = float(W[0].sum())
A = ROW_SUM * np.eye(N) - W
EVALS, EVECS = np.linalg.eigh(A)


def spectral_function(values: np.ndarray) -> np.ndarray:
    return (EVECS * values) @ EVECS.T


def markov(t: float) -> np.ndarray:
    return spectral_function(np.exp(-t * EVALS)).real


def unitary(t: float) -> np.ndarray:
    return spectral_function(np.exp(-1j * t * EVALS))


def distribution_diffusive(t: float) -> np.ndarray:
    p = markov(t)[:, 0]
    p[np.abs(p) < 1e-15] = 0.0
    return p / p.sum()


def distribution_unitary(t: float) -> np.ndarray:
    q = np.abs(unitary(t)[:, 0]) ** 2
    return q / q.sum()


def entropy(p: np.ndarray) -> float:
    x = np.asarray(p, dtype=float)
    x = x[x > 0.0]
    return float(-(x * np.log(x)).sum())


def tv(p: np.ndarray, q: np.ndarray) -> float:
    return float(0.5 * np.abs(np.asarray(p) - np.asarray(q)).sum())


def chernoff(p: np.ndarray, q: np.ndarray) -> float:
    p = np.maximum(np.asarray(p, dtype=float), 1e-300)
    q = np.maximum(np.asarray(q, dtype=float), 1e-300)

    def objective(alpha: float) -> float:
        return float(np.sum(np.exp(alpha * np.log(p) + (1.0 - alpha) * np.log(q))))

    lo, hi = 0.0, 1.0
    for _ in range(70):
        a = (2.0 * lo + hi) / 3.0
        b = (lo + 2.0 * hi) / 3.0
        if objective(a) <= objective(b):
            hi = b
        else:
            lo = a
    return float(-math.log(objective((lo + hi) / 2.0)))


def population_matrix(t: float) -> np.ndarray:
    return np.abs(unitary(t)) ** 2


def joint_record(transition: np.ndarray) -> np.ndarray:
    first = transition[:, 0]
    # joint[x1, x2] = Pr(x1) Pr(x2 | x1)
    return first[:, None] * transition.T


def optimize_grid(
    objective: Callable[[float], float], start: float, stop: float, step: float
) -> tuple[float, float]:
    best_t, best_value = start, -math.inf
    count = int(round((stop - start) / step))
    for idx in range(count + 1):
        t = start + idx * step
        value = objective(t)
        if value > best_value:
            best_t, best_value = t, value
    return float(best_t), float(best_value)


def apparatus_distribution(base: np.ndarray, fs: tuple[int, int], fa: tuple[int, int]) -> np.ndarray:
    rs, _ls = fs
    ra, la = fa
    return np.array([base[(ra + la * delta - rs) % N] for delta in range(N)])


def transform_frame(g: tuple[int, int], f: tuple[int, int]) -> tuple[int, int]:
    a, eps = g
    r, lam = f
    return ((a + eps * r) % N, eps * lam)


def relational_covariance_residual(base: np.ndarray) -> float:
    frames = [(r, lam) for r in range(N) for lam in (-1, 1)]
    group = [(a, eps) for a in range(N) for eps in (-1, 1)]
    maximum = 0.0
    for fs in frames:
        for fa in frames:
            reference = apparatus_distribution(base, fs, fa)
            for g in group:
                transformed = apparatus_distribution(
                    base, transform_frame(g, fs), transform_frame(g, fa)
                )
                maximum = max(maximum, float(np.max(np.abs(reference - transformed))))
    return maximum


def odd_probe() -> dict[str, float]:
    shift = np.zeros((N, N), dtype=complex)
    for i in range(N):
        shift[(i + 1) % N, i] = 1.0
    current = (shift - shift.conj().T) / (2.0j)
    h = 0.1
    hamiltonian = W - h * current
    vals, vecs = np.linalg.eigh(hamiltonian)
    u = (vecs * np.exp(-2.0j * vals)) @ vecs.conj().T
    probs = np.abs(u[:, 0]) ** 2
    p_plus = float(probs[1])
    p_minus = float(probs[-1])
    # Fix the apparatus name so that plus is the favored direction.
    if p_plus < p_minus:
        p_plus, p_minus = p_minus, p_plus
    accuracy = p_plus / (p_plus + p_minus)

    result: dict[str, float] = {
        "p_plus": p_plus,
        "p_minus": p_minus,
        "conditional_accuracy": accuracy,
    }
    for reliability in (0.5, 0.6, 0.8, 0.95, 1.0):
        effective = reliability * accuracy + (1.0 - reliability) * (1.0 - accuracy)
        if effective in (0.0, 1.0):
            mutual = 1.0
        else:
            h2 = -effective * math.log2(effective) - (1.0 - effective) * math.log2(
                1.0 - effective
            )
            mutual = 1.0 - h2
        result[f"mutual_information_r={reliability:.2f}"] = mutual
    return result


def main() -> None:
    assert np.max(np.abs(W.sum(axis=1) - ROW_SUM)) < 1e-14
    assert float(EVALS.min()) > -1e-12
    assert np.min(W + np.diag(np.full(N, math.inf))) > 0.0

    basic: dict[str, object] = {
        "row_sum": ROW_SUM,
        "rates": [float(x) for x in EVALS],
        "dirichlet_min_eigenvalue": float(EVALS.min()),
        "variance_short_time": float(np.sum(W[:, 0] ** 2)),
        "snapshots": {},
    }

    for t in (0.1, 1.0, 5.0):
        q = distribution_unitary(t)
        p = distribution_diffusive(t)
        basic["snapshots"][str(t)] = {
            "unitary_return": float(q[0]),
            "unitary_entropy": entropy(q),
            "markov_return": float(p[0]),
            "markov_entropy": entropy(p),
            "total_variation": tv(q, p),
        }

    tv_time, tv_best = optimize_grid(
        lambda t: tv(distribution_unitary(t), distribution_diffusive(t)),
        0.005,
        10.0,
        0.005,
    )
    ch_time, ch_best = optimize_grid(
        lambda t: chernoff(distribution_unitary(t), distribution_diffusive(t)),
        0.005,
        10.0,
        0.005,
    )

    ck: dict[str, object] = {}
    for t in (0.05, 0.1, 0.5, 1.0, 2.0, 3.0):
        pt = markov(t)
        mt = population_matrix(t)
        ck[str(t)] = {
            "markov_frobenius": float(np.linalg.norm(markov(2.0 * t) - pt @ pt)),
            "unitary_population_frobenius": float(
                np.linalg.norm(population_matrix(2.0 * t) - mt @ mt)
            ),
            "unitary_population_tv_from_origin": tv(
                population_matrix(2.0 * t)[:, 0], (mt @ mt)[:, 0]
            ),
        }

    record: dict[str, object] = {}
    for t in (0.25, 0.5, 1.0, 2.0, 3.0, 5.0):
        jp = joint_record(markov(t))
        jq = joint_record(population_matrix(t))
        record[str(t)] = {
            "joint_tv": tv(jp.ravel(), jq.ravel()),
            "terminal_tv": tv(jp.sum(axis=0), jq.sum(axis=0)),
            "bayes_error": 0.5 * (1.0 - tv(jp.ravel(), jq.ravel())),
        }

    record_time, record_best = optimize_grid(
        lambda t: tv(
            joint_record(markov(t)).ravel(),
            joint_record(population_matrix(t)).ravel(),
        ),
        0.005,
        5.0,
        0.005,
    )

    zeno: dict[str, object] = {}
    for total, repetitions in ((2.0, 10), (2.0, 200), (2.0, 1000), (5.0, 10), (5.0, 1000)):
        delta = total / repetitions
        m = population_matrix(delta)
        quantum_survival = float(distribution_unitary(delta)[0] ** repetitions)
        classical_checked = float(distribution_diffusive(delta)[0] ** repetitions)
        dephased = np.linalg.matrix_power(m, repetitions)[:, 0]
        zeno[f"T={total},n={repetitions}"] = {
            "quantum_survival": quantum_survival,
            "uninterrupted_unitary": float(distribution_unitary(total)[0]),
            "classical_no_jump_proxy": classical_checked,
            "analytic_classical_no_jump": math.exp(-ROW_SUM * total),
            "dephased_origin": float(dephased[0]),
            "dephased_tv_to_diffusion": tv(dephased, distribution_diffusive(total)),
        }

    times = np.arange(0.0, 20.0 + 0.005, 0.005)
    ent_u = np.array([entropy(distribution_unitary(float(t))) for t in times])
    ent_d = np.array([entropy(distribution_diffusive(float(t))) for t in times])
    entropy_scan = {
        "unitary_decreasing_steps": int(np.sum(np.diff(ent_u) < -1e-12)),
        "diffusive_decreasing_steps": int(np.sum(np.diff(ent_d) < -1e-12)),
        "markov_entropy_t=20": float(ent_d[-1]),
        "log_12": math.log(12.0),
    }

    cycle_scaling: dict[str, object] = {}
    for n in (12, 24, 48, 96, 192):
        gap = 4.0 * math.sin(math.pi / n) ** 2
        cycle_scaling[str(n)] = {
            "gap": gap,
            "diffusive_time": 1.0 / gap,
            "wave_time": 1.0 / math.sqrt(gap),
        }

    base_u = distribution_unitary(2.0)
    base_d = distribution_diffusive(2.0)
    relational = {
        "unitary_covariance_residual": relational_covariance_residual(base_u),
        "diffusive_covariance_residual": relational_covariance_residual(base_d),
        "odd_probe": odd_probe(),
    }

    result = {
        "dual_generator": basic,
        "optimal_single_time_tv": {"time": tv_time, "value": tv_best},
        "optimal_single_time_chernoff": {"time": ch_time, "value": ch_best},
        "chapman_kolmogorov": ck,
        "two_step_records": record,
        "optimal_record_tv": {"time": record_time, "value": record_best},
        "zeno_and_dephasing": zeno,
        "entropy_scan": entropy_scan,
        "cycle_scaling": cycle_scaling,
        "relational_apparatus": relational,
    }

    assert max(item["markov_frobenius"] for item in ck.values()) < 1e-12
    assert ck["0.5"]["unitary_population_frobenius"] > 0.5
    assert entropy_scan["diffusive_decreasing_steps"] == 0
    assert relational["unitary_covariance_residual"] < 1e-12
    assert relational["diffusive_covariance_residual"] < 1e-12

    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

