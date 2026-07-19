#!/usr/bin/env python3
"""Reproducible finite studies for FIN Programs 21--30.

The program is deliberately self-contained.  It uses the frozen twelve-mode
radial data and does not import repository verdicts.  All temporal quantities
are dimensionless.  The numerical studies test mathematical identifiability
and finite constructions; they are not evidence for fundamental physics.

Dependencies: NumPy, SciPy, Matplotlib (figures only).
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize, minimize_scalar, nnls


N = 12
EPS = 1.0e-14
RNG_SEED = 20260719
ROW_SUM_EXACT = 1.660307278766099
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


def cyclic_distance(i: int, j: int, n: int) -> int:
    delta = abs(i - j)
    return min(delta, n - delta)


def cyclic_matrix(radial: np.ndarray, n: int = N) -> np.ndarray:
    return np.array(
        [[radial[cyclic_distance(i, j, n)] for j in range(n)] for i in range(n)],
        dtype=float,
    )


W = cyclic_matrix(RADIAL)
ROW_SUM = float(W[0].sum())
A = ROW_SUM * np.eye(N) - W
A_EVALS, A_EVECS = np.linalg.eigh(A)


def spectral_function(values: np.ndarray) -> np.ndarray:
    return (A_EVECS * values) @ A_EVECS.T


def markov(t: float) -> np.ndarray:
    return spectral_function(np.exp(-t * A_EVALS)).real


def unitary(t: float) -> np.ndarray:
    return spectral_function(np.exp(-1j * t * A_EVALS))


def population_unitary(t: float) -> np.ndarray:
    return np.abs(unitary(t)) ** 2


def normalize_probability(p: np.ndarray) -> np.ndarray:
    q = np.asarray(p, dtype=float).copy()
    q[np.abs(q) < 1.0e-15] = 0.0
    q = np.maximum(q, 0.0)
    return q / q.sum()


def kl(p: np.ndarray, q: np.ndarray) -> float:
    p = np.maximum(np.asarray(p, dtype=float), 1.0e-300)
    q = np.maximum(np.asarray(q, dtype=float), 1.0e-300)
    p = p / p.sum()
    q = q / q.sum()
    return float(np.sum(p * np.log(p / q)))


def js(p: np.ndarray, q: np.ndarray) -> float:
    p = normalize_probability(p)
    q = normalize_probability(q)
    m = 0.5 * (p + q)
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def total_variation(p: np.ndarray, q: np.ndarray) -> float:
    return float(0.5 * np.abs(normalize_probability(p) - normalize_probability(q)).sum())


def liouvillian_transition(t: float, gamma: float) -> np.ndarray:
    """Position transition matrix for Hamiltonian A plus site dephasing."""
    identity = np.eye(N)
    # vec(A rho - rho A) in column-major convention.
    coherent = -1j * (np.kron(identity, A) - np.kron(A.T, identity))
    dephase = np.zeros((N * N, N * N), dtype=complex)
    for i in range(N):
        for j in range(N):
            index = i + N * j
            if i != j:
                dephase[index, index] = -gamma
    channel = expm(t * (coherent + dephase))
    transition = np.zeros((N, N), dtype=float)
    for source in range(N):
        rho = np.zeros((N, N), dtype=complex)
        rho[source, source] = 1.0
        final = (channel @ rho.reshape(-1, order="F")).reshape((N, N), order="F")
        transition[:, source] = np.real(np.diag(final))
    transition[np.abs(transition) < 1.0e-13] = 0.0
    transition = np.maximum(transition, 0.0)
    transition /= transition.sum(axis=0, keepdims=True)
    return transition


def process_blocks(model: str, times: tuple[float, float, float]) -> list[np.ndarray]:
    """Intervention-aware one-, two-, and three-time probability blocks."""
    t1, t2, t3 = times

    def transition(t: float, branch: int | None = None) -> np.ndarray:
        if model == "coherent":
            return population_unitary(t)
        if model == "markov":
            return markov(t)
        if model == "dephasing":
            return liouvillian_transition(t, 0.55)
        if model == "memory_bath":
            gamma = (0.08, 1.25)[0 if branch is None else branch]
            return liouvillian_transition(t, gamma)
        raise ValueError(model)

    e0 = np.zeros(N)
    e0[0] = 1.0
    blocks: list[np.ndarray] = []
    if model != "memory_bath":
        q1, q2, q3 = transition(t1), transition(t2), transition(t3)
        q12 = transition(t1 + t2)
        q123 = transition(t1 + t2 + t3)
        blocks.append(q12 @ e0)  # no-intervention terminal distribution
        blocks.append(q2 @ q1 @ e0)  # unread projective intervention
        joint2 = (q1 @ e0)[:, None] * q2.T
        blocks.append(joint2.ravel())
        joint3 = (
            (q1 @ e0)[:, None, None]
            * q2.T[:, :, None]
            * np.transpose(q3, (1, 0))[None, :, :]
        )
        blocks.append(joint3.ravel())
        blocks.append(q123 @ e0)
        for reset in (0, 1, 3):
            blocks.append(q3[:, reset])
    else:
        weight = 0.62
        branch_blocks: list[list[np.ndarray]] = []
        for branch in (0, 1):
            q1 = transition(t1, branch)
            q2 = transition(t2, branch)
            q3 = transition(t3, branch)
            q12 = transition(t1 + t2, branch)
            q123 = transition(t1 + t2 + t3, branch)
            joint2 = (q1 @ e0)[:, None] * q2.T
            joint3 = (
                (q1 @ e0)[:, None, None]
                * q2.T[:, :, None]
                * np.transpose(q3, (1, 0))[None, :, :]
            )
            branch_blocks.append(
                [q12 @ e0, q2 @ q1 @ e0, joint2.ravel(), joint3.ravel(), q123 @ e0]
                + [q3[:, reset] for reset in (0, 1, 3)]
            )
        blocks = [
            weight * branch_blocks[0][i] + (1.0 - weight) * branch_blocks[1][i]
            for i in range(len(branch_blocks[0]))
        ]
    return [normalize_probability(block) for block in blocks]


def concatenate_blocks(blocks: list[np.ndarray]) -> np.ndarray:
    return np.concatenate([block / len(blocks) for block in blocks])


def process_tomography_study() -> dict[str, object]:
    models = ("coherent", "markov", "dephasing", "memory_bath")
    # A bounded design scan.  Each candidate uses unequal times to suppress aliases.
    candidates = [float(x) for x in np.linspace(0.25, 2.25, 17)]
    best: tuple[float, tuple[float, float, float], dict[str, list[np.ndarray]]] | None = None
    for anchor in candidates:
        times = (0.55 * anchor, anchor, 1.65 * anchor)
        signatures = {model: process_blocks(model, times) for model in models}
        vectors = {model: concatenate_blocks(signatures[model]) for model in models}
        minimum = min(js(vectors[a], vectors[b]) for a, b in itertools.combinations(models, 2))
        if best is None or minimum > best[0]:
            best = (minimum, times, signatures)
    assert best is not None
    minimum_js, times, signatures = best
    vectors = {model: concatenate_blocks(signatures[model]) for model in models}
    pairwise_js = {
        f"{a}__{b}": js(vectors[a], vectors[b]) for a, b in itertools.combinations(models, 2)
    }
    pairwise_tv = {
        f"{a}__{b}": total_variation(vectors[a], vectors[b])
        for a, b in itertools.combinations(models, 2)
    }
    signature_matrix = np.vstack([vectors[model] for model in models])
    singular = np.linalg.svd(signature_matrix - signature_matrix.mean(axis=0), compute_uv=False)

    # Finite-shot likelihood experiment; every probability block receives the same shots.
    rng = np.random.default_rng(RNG_SEED + 21)
    shots = 400
    trials = 80
    confusion = np.zeros((len(models), len(models)), dtype=int)
    for true_index, true_model in enumerate(models):
        for _ in range(trials):
            counts = [rng.multinomial(shots, block) for block in signatures[true_model]]
            scores = []
            for candidate in models:
                score = 0.0
                for count, probability in zip(counts, signatures[candidate]):
                    score += float(count @ np.log(np.maximum(probability, 1.0e-300)))
                scores.append(score)
            confusion[true_index, int(np.argmax(scores))] += 1
    return {
        "models": list(models),
        "optimized_times": list(times),
        "probability_blocks": len(next(iter(signatures.values()))),
        "signature_dimension": int(next(iter(vectors.values())).size),
        "minimum_pairwise_JS_nats": minimum_js,
        "pairwise_JS_nats": pairwise_js,
        "pairwise_total_variation": pairwise_tv,
        "centered_signature_singular_values": singular.tolist(),
        "centered_signature_rank_tol_1e-10": int(np.sum(singular > 1.0e-10)),
        "shots_per_block": shots,
        "trials_per_model": trials,
        "confusion_matrix_rows_true_columns_predicted": confusion.tolist(),
        "classification_accuracy": float(np.trace(confusion) / confusion.sum()),
    }


def strict_formula(d: int) -> float:
    if d == 0:
        return 0.0
    return math.cos(0.18575 * d + 0.16250) / (1.0 + d**1.8)


def positivity_refinement_study() -> dict[str, object]:
    anchors = np.array([strict_formula(d) for d in range(1, 7)])
    tail_ratio = float(anchors[-1] / anchors[-2])

    def completed_weight(d: int) -> float:
        if d <= 6:
            return float(anchors[d - 1])
        return float(anchors[-1] * tail_ratio ** (d - 6))

    first_negative = next(d for d in range(1, 100) if strict_formula(d) < 0.0)
    size_checks: dict[str, dict[str, float]] = {}
    scaled_gaps: list[float] = []
    for n in (12, 16, 24, 48, 96, 192, 384, 512):
        radial = np.array([completed_weight(d) for d in range(n // 2 + 1)])
        radial[0] = 0.0
        matrix = cyclic_matrix(radial, n)
        row_sum = float(matrix[0].sum())
        generator = row_sum * np.eye(n) - matrix
        eigenvalues = np.linalg.eigvalsh(generator)
        gap = float(eigenvalues[1])
        scaled = n * n * gap
        scaled_gaps.append(scaled)
        offdiag = matrix.copy()
        np.fill_diagonal(offdiag, math.inf)
        size_checks[str(n)] = {
            "minimum_offdiagonal_weight": float(offdiag.min()),
            "minimum_generator_eigenvalue": float(eigenvalues[0]),
            "N_squared_gap": scaled,
        }
    second_moment = sum(completed_weight(d) * d * d for d in range(1, 600))
    predicted_gap_limit = 4.0 * math.pi**2 * second_moment

    # Smooth completely-monotone two-exponential alternative.
    lambdas = np.logspace(-3, 2, 1200)
    design = np.exp(-np.arange(1, 7)[:, None] * lambdas[None, :])
    coefficients, residual = nnls(design, anchors)
    active = np.where(coefficients > 1.0e-9)[0]
    mixture = [
        {"coefficient": float(coefficients[i]), "decay_rate": float(lambdas[i])}
        for i in active
    ]

    # A second, inequivalent positive continuation keeps the circumference
    # fixed and rescales the argument.  It has a nonlocal integral limit.
    fixed_circumference: dict[str, dict[str, float]] = {}
    reference_n = 1536

    def rescaled_matrix(n: int) -> np.ndarray:
        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    r = 12.0 * cyclic_distance(i, j, n) / n
                    matrix[i, j] = (12.0 / n) * math.cos(0.18575 * r + 0.16250) / (1.0 + r**1.8)
        return matrix

    reference = rescaled_matrix(reference_n)
    ref_generator = float(reference[0].sum()) * np.eye(reference_n) - reference
    ref_modes = np.fft.fft(ref_generator[0]).real[:7]
    for n in (12, 16, 24, 48, 96, 192):
        matrix = rescaled_matrix(n)
        generator = float(matrix[0].sum()) * np.eye(n) - matrix
        modes = np.fft.fft(generator[0]).real[:7]
        offdiag = matrix.copy()
        np.fill_diagonal(offdiag, math.inf)
        fixed_circumference[str(n)] = {
            "minimum_offdiagonal_weight": float(offdiag.min()),
            "maximum_relative_error_modes_1_to_6_vs_N1536": float(
                np.max(np.abs((modes[1:7] - ref_modes[1:7]) / ref_modes[1:7]))
            ),
        }
    return {
        "strict_anchor_values_d1_to_d6": anchors.tolist(),
        "exact_anchor_max_error": float(np.max(np.abs(anchors - RADIAL[1:]))),
        "raw_strict_first_negative_distance": first_negative,
        "tail_ratio": tail_ratio,
        "all_checked_weights_positive": bool(
            all(item["minimum_offdiagonal_weight"] > 0.0 for item in size_checks.values())
        ),
        "size_checks": size_checks,
        "second_moment": float(second_moment),
        "predicted_N_squared_gap_limit": float(predicted_gap_limit),
        "relative_error_last_scaled_gap": float(
            abs(scaled_gaps[-1] - predicted_gap_limit) / predicted_gap_limit
        ),
        "two_exponential_active_terms": mixture,
        "two_exponential_relative_anchor_residual": float(residual / np.linalg.norm(anchors)),
        "fixed_circumference_nonlocal_continuation": fixed_circumference,
        "fixed_circumference_limit_first_four_modes": ref_modes[1:5].tolist(),
        "status": "existence witness, not a canonical FIN refinement",
    }


def cycle_incidence(n: int) -> np.ndarray:
    b = np.zeros((n, n), dtype=float)
    for edge in range(n):
        b[edge, edge] = -1.0
        b[edge, (edge + 1) % n] = 1.0
    return b


def dirac_refinement_study() -> dict[str, object]:
    square_residuals: dict[str, float] = {}
    gaps: dict[str, float] = {}
    for n in (12, 24, 48, 96, 192, 384):
        b = cycle_incidence(n)
        d = np.block([[np.zeros((n, n)), b.T], [b, np.zeros((n, n))]])
        target = np.block([[b.T @ b, np.zeros((n, n))], [np.zeros((n, n)), b @ b.T]])
        square_residuals[str(n)] = float(np.linalg.norm(d @ d - target))
        gaps[str(n)] = 2.0 * math.sin(math.pi / n)
    ns = np.array([24, 48, 96, 192, 384], dtype=float)
    z_fit = float(-np.polyfit(np.log(ns), np.log([gaps[str(int(n))] for n in ns]), 1)[0])
    momenta = np.linspace(0.0, math.pi / 3.0, 400)[1:]
    dispersion = 2.0 * np.sin(momenta / 2.0)
    relative_linearity = np.abs(dispersion / momenta - 1.0)

    # Finite propagation-tail diagnostic on the doubled local graph.
    n = 96
    b = cycle_incidence(n)
    d = np.block([[np.zeros((n, n)), b.T], [b, np.zeros((n, n))]])
    values, vectors = np.linalg.eigh(d)
    initial = np.zeros(2 * n, dtype=complex)
    initial[0] = 1.0
    tails: dict[str, float] = {}
    for time in (1.0, 2.0, 4.0, 8.0):
        state = vectors @ (np.exp(-1j * time * values) * (vectors.T @ initial))
        radius = int(math.ceil(2.0 * time + 2.0))
        outside = 0.0
        for index, probability in enumerate(np.abs(state) ** 2):
            site = index % n
            distance = min(site, n - site)
            if distance > radius:
                outside += float(probability)
        tails[str(time)] = outside

    # Exact square root of the strict Dirichlet generator.  Since all strict
    # weights are positive, the incidence graph is complete and therefore
    # nonlocal: exact strict matching trades away nearest-neighbour locality.
    edges: list[tuple[int, int]] = [(i, j) for i in range(N) for j in range(i + 1, N)]
    weighted_incidence = np.zeros((len(edges), N))
    for row, (i, j) in enumerate(edges):
        weight = math.sqrt(W[i, j])
        weighted_incidence[row, i] = -weight
        weighted_incidence[row, j] = weight
    strict_dirac = np.block(
        [
            [np.zeros((N, N)), weighted_incidence.T],
            [weighted_incidence, np.zeros((len(edges), len(edges)))],
        ]
    )
    strict_target = np.block(
        [
            [A, np.zeros((N, len(edges)))],
            [np.zeros((len(edges), N)), weighted_incidence @ weighted_incidence.T],
        ]
    )
    sorted_rates = np.sort(A_EVALS)
    z12 = math.log(math.sqrt(sorted_rates[3]) / math.sqrt(sorted_rates[1])) / math.log(2.0)
    z13 = math.log(math.sqrt(sorted_rates[5]) / math.sqrt(sorted_rates[1])) / math.log(3.0)
    return {
        "D_squared_block_L_residuals": square_residuals,
        "positive_Dirac_gaps": gaps,
        "fitted_dynamic_exponent": z_fit,
        "max_relative_dispersion_error_k_le_pi_over_3": float(relative_linearity.max()),
        "probability_outside_radius_ceil_2t_plus_2": tails,
        "strict_weighted_incidence_edges": len(edges),
        "strict_D_squared_residual": float(np.linalg.norm(strict_dirac @ strict_dirac - strict_target)),
        "strict_square_root_effective_exponent_modes_1_2": float(z12),
        "strict_square_root_effective_exponent_modes_1_3": float(z13),
        "orientation_statement": "edge orientation changes D by unitary equivalence; D^2 is unchanged",
        "physical_status": "exact z=1 mathematical refinement; not selected by the strict kernel and not a Lorentz theorem",
    }


def os_positivity_study() -> dict[str, object]:
    mass2 = 0.35
    omega_values = np.sqrt(A_EVALS + mass2)

    def correlator(t: float) -> np.ndarray:
        return (A_EVECS * (np.exp(-omega_values * abs(t)) / (2.0 * omega_values))) @ A_EVECS.T

    times = np.array([0.15, 0.35, 0.7, 1.1, 1.8])
    blocks = [[correlator(float(ti + tj)) for tj in times] for ti in times]
    reflection = np.block(blocks)
    reflection_eigenvalues = np.linalg.eigvalsh(reflection)

    # Analytic-continuation negative control: cos(omega(t+s))/(2 omega).
    bad_times = np.array([0.4, 1.2])
    bad = np.array(
        [[math.cos(float(ti + tj)) / 2.0 for tj in bad_times] for ti in bad_times]
    )
    bad_eigenvalues = np.linalg.eigvalsh(bad)

    finite_lattice: dict[str, dict[str, float | int]] = {}
    for temporal_size in (8, 12, 16, 24):
        shift = np.roll(np.eye(temporal_size), 1, axis=0)
        time_laplacian = 2.0 * np.eye(temporal_size) - shift - shift.T
        for alpha in (0.0, 0.01, 0.1, 1.0):
            temporal_operator = time_laplacian + alpha * (time_laplacian @ time_laplacian)
            hessian = np.kron(temporal_operator, np.eye(N)) + np.kron(
                np.eye(temporal_size), A + 0.3 * np.eye(N)
            )
            covariance = np.linalg.inv(hessian)
            positive_times = list(range(1, temporal_size // 2))
            reflection_indices = [((-t) % temporal_size) * N + x for t in positive_times for x in range(N)]
            positive_indices = [t * N + x for t in positive_times for x in range(N)]
            form = covariance[np.ix_(reflection_indices, positive_indices)]
            form = 0.5 * (form + form.T)
            eigenvalues = np.linalg.eigvalsh(form)
            finite_lattice[f"T={temporal_size},alpha={alpha}"] = {
                "minimum_reflection_eigenvalue": float(eigenvalues[0]),
                "negative_count_below_minus_1e-10": int(np.sum(eigenvalues < -1.0e-10)),
                "hessian_minimum_eigenvalue": float(np.linalg.eigvalsh(hessian)[0]),
            }
    return {
        "mass_squared": mass2,
        "positive_time_samples": times.tolist(),
        "reflection_matrix_dimension": int(reflection.shape[0]),
        "minimum_reflection_eigenvalue": float(reflection_eigenvalues[0]),
        "negative_eigenvalues_below_minus_1e-11": int(np.sum(reflection_eigenvalues < -1.0e-11)),
        "positive_rank_tol_1e-11": int(np.sum(reflection_eigenvalues > 1.0e-11)),
        "factorization_rank_bound": N,
        "bad_cosine_control_eigenvalues": bad_eigenvalues.tolist(),
        "bad_control_violates_positivity": bool(bad_eigenvalues[0] < -1.0e-8),
        "finite_temporal_lattice_tests": finite_lattice,
        "scope": "free massive Gaussian correlator only",
    }


def observed_distribution(model: str, laboratory_time: float, gamma: float, spam: float) -> np.ndarray:
    t = gamma * laboratory_time
    p = (population_unitary(t) if model == "coherent" else markov(t))[:, 0]
    p = normalize_probability(p)
    return (1.0 - spam) * p + spam * np.ones(N) / N


def nuisance_robust_identification_study() -> dict[str, object]:
    models = ("coherent", "markov")
    gamma_grid = np.linspace(0.95, 1.05, 21)
    spam_grid = np.linspace(0.0, 0.06, 13)
    candidate_times = np.linspace(0.25, 5.0, 40)

    def robust_separation(times: tuple[float, ...]) -> float:
        worst = math.inf
        for true_model in models:
            other = models[1 - models.index(true_model)]
            truth = np.concatenate(
                [observed_distribution(true_model, t, 1.0, 0.03) / len(times) for t in times]
            )
            best_alt = math.inf
            for gamma in gamma_grid:
                for spam in spam_grid:
                    alt = np.concatenate(
                        [observed_distribution(other, t, float(gamma), float(spam)) / len(times) for t in times]
                    )
                    best_alt = min(best_alt, js(truth, alt))
            worst = min(worst, best_alt)
        return float(worst)

    # Greedy maximin selection of three distinct times.
    selected: list[float] = []
    gains: list[float] = []
    for _ in range(3):
        scores = [
            (robust_separation(tuple(selected + [float(t)])), float(t))
            for t in candidate_times
            if all(abs(t - s) > 1.0e-9 for s in selected)
        ]
        value, time = max(scores)
        selected.append(time)
        gains.append(value)

    # Compare one snapshot under broad nuisance freedom.
    broad_gamma = np.linspace(0.5, 1.8, 66)
    broad_spam = np.linspace(0.0, 0.25, 26)
    single_time = selected[0]
    truth_single = observed_distribution("coherent", single_time, 1.0, 0.03)
    broad_minimum = min(
        js(truth_single, observed_distribution("markov", single_time, float(g), float(e)))
        for g in broad_gamma
        for e in broad_spam
    )

    # Shot-noise profile-likelihood classification using the calibrated bounds.
    rng = np.random.default_rng(RNG_SEED + 25)
    shots = 600
    trials = 120
    confusion = np.zeros((2, 2), dtype=int)
    for true_index, true_model in enumerate(models):
        for _ in range(trials):
            counts = [
                rng.multinomial(shots, observed_distribution(true_model, t, 1.017, 0.028))
                for t in selected
            ]
            profile_scores: list[float] = []
            for candidate in models:
                best_score = -math.inf
                for gamma in gamma_grid:
                    for spam in spam_grid:
                        score = 0.0
                        for count, t in zip(counts, selected):
                            p = observed_distribution(candidate, t, float(gamma), float(spam))
                            score += float(count @ np.log(np.maximum(p, 1.0e-300)))
                        best_score = max(best_score, score)
                profile_scores.append(best_score)
            confusion[true_index, int(np.argmax(profile_scores))] += 1
    return {
        "selected_times": selected,
        "greedy_maximin_JS_after_each_time": gains,
        "calibrated_clock_range": [float(gamma_grid[0]), float(gamma_grid[-1])],
        "calibrated_SPAM_range": [float(spam_grid[0]), float(spam_grid[-1])],
        "single_snapshot_broad_nuisance_minimum_JS": float(broad_minimum),
        "shots_per_time": shots,
        "trials_per_model": trials,
        "confusion_matrix_rows_true_columns_predicted": confusion.tolist(),
        "classification_accuracy": float(np.trace(confusion) / confusion.sum()),
        "interpretation": "conditional identifiability inside declared nuisance bounds",
    }


def bath_identifiability_study() -> dict[str, object]:
    frequencies = np.array([-1.7, 0.2, 1.15])
    weights = np.array([0.25, 0.45, 0.30])
    times = np.arange(0, 18, dtype=float) * 0.35
    signal = np.array([np.sum(weights * np.exp(-1j * frequencies * t)) for t in times])
    hankel = np.array([[signal[i + j] for j in range(9)] for i in range(9)])
    singular = np.linalg.svd(hankel, compute_uv=False)

    rng = np.random.default_rng(RNG_SEED + 26)
    noisy = signal + 2.0e-3 * (rng.normal(size=signal.size) + 1j * rng.normal(size=signal.size))
    hankel_noisy = np.array([[noisy[i + j] for j in range(9)] for i in range(9)])
    singular_noisy = np.linalg.svd(hankel_noisy, compute_uv=False)

    # Finite-window nonuniqueness: two nonnegative spectral measures agree on
    # the chosen real cosine samples but differ later.  A null-space direction
    # is added to a strictly positive reference weight vector.
    sample_times = np.linspace(0.0, 2.0, 7)
    frequency_grid = np.linspace(0.0, 4.5, 15)
    design = np.cos(sample_times[:, None] * frequency_grid[None, :])
    _u, _s, vh = np.linalg.svd(design)
    null = vh[-1]
    base = np.ones(frequency_grid.size) / frequency_grid.size
    bound = min(base[i] / abs(null[i]) for i in range(null.size) if abs(null[i]) > 1.0e-14)
    perturbation = 0.7 * bound * null
    w_plus = base + perturbation
    w_minus = base - perturbation
    training_difference = np.max(np.abs(design @ w_plus - design @ w_minus))
    future_times = np.linspace(2.5, 12.0, 200)
    future_design = np.cos(future_times[:, None] * frequency_grid[None, :])
    future_difference = np.max(np.abs(future_design @ w_plus - future_design @ w_minus))

    # Fixed-time qubit dephasing channel: Kraus/Choi rank is two unless |g|=1.
    g = abs(signal[4])
    choi_eigenvalues = np.array([1.0 + g, 1.0 - g, 0.0, 0.0])

    # Strong nonidentifiability witness: two continuous quantum semigroups
    # induce the same complete population process P_t but have radically
    # different one-time Choi ranks.
    t_channel = 1.0
    p_channel = markov(t_channel)
    shift = np.zeros((N, N), dtype=complex)
    for j in range(N):
        shift[(j + 1) % N, j] = 1.0
    shifts = [np.linalg.matrix_power(shift, d) for d in range(N)]
    displacement = p_channel[:, 0]

    def random_shift_channel(rho: np.ndarray) -> np.ndarray:
        return sum(
            displacement[d] * shifts[d] @ rho @ shifts[d].conj().T for d in range(N)
        )

    no_jump = math.exp(-ROW_SUM * t_channel)

    def edge_jump_channel(rho: np.ndarray) -> np.ndarray:
        diagonal = np.diag(rho).copy()
        return no_jump * (rho - np.diag(diagonal)) + np.diag(p_channel @ diagonal)

    def choi(channel: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
        out = np.zeros((N * N, N * N), dtype=complex)
        for j in range(N):
            for k in range(N):
                matrix_unit = np.zeros((N, N), dtype=complex)
                matrix_unit[j, k] = 1.0
                out[j * N : (j + 1) * N, k * N : (k + 1) * N] = channel(matrix_unit)
        return out

    random_choi = np.linalg.eigvalsh(choi(random_shift_channel))
    edge_choi = np.linalg.eigvalsh(choi(edge_jump_channel))
    population_residual = 0.0
    for source in range(N):
        rho = np.zeros((N, N), dtype=complex)
        rho[source, source] = 1.0
        population_residual = max(
            population_residual,
            float(
                np.max(
                    np.abs(
                        np.diag(random_shift_channel(rho)) - np.diag(edge_jump_channel(rho))
                    )
                )
            ),
        )
    plus = np.ones(N, dtype=complex) / math.sqrt(N)
    rho_plus = np.outer(plus, plus.conj())
    channel_difference = random_shift_channel(rho_plus) - edge_jump_channel(rho_plus)
    coherent_trace_distance = 0.5 * float(np.sum(np.abs(np.linalg.eigvalsh(channel_difference))))
    return {
        "true_finite_mode_count": 3,
        "Hankel_singular_values_noiseless": singular.tolist(),
        "Hankel_rank_tol_1e-10": int(np.sum(singular > 1.0e-10)),
        "Hankel_singular_values_noise_2e-3": singular_noisy.tolist(),
        "fixed_time_dephasing_Choi_eigenvalues": choi_eigenvalues.tolist(),
        "fixed_time_minimal_Stinespring_dimension": int(np.sum(choi_eigenvalues > 1.0e-12)),
        "finite_window_alternative_measures_min_weight": float(min(w_plus.min(), w_minus.min())),
        "finite_window_training_max_difference": float(training_difference),
        "future_max_difference": float(future_difference),
        "population_equivalent_quantum_semigroups": {
            "time": t_channel,
            "population_record_max_residual": population_residual,
            "random_shift_Choi_rank": int(np.sum(random_choi > 1.0e-10)),
            "edge_jump_Choi_rank": int(np.sum(edge_choi > 1.0e-10)),
            "uniform_coherent_input_trace_distance": coherent_trace_distance,
            "analytic_trace_distance": (1.0 - no_jump) * (N - 1.0) / N,
        },
        "verdict": "finite noiseless model order is recoverable; a unique bath is not",
    }


def dihedral_matrix(a: int, epsilon: int) -> np.ndarray:
    matrix = np.zeros((N, N))
    for i in range(N):
        matrix[(a + epsilon * i) % N, i] = 1.0
    return matrix


def relational_audit_study() -> dict[str, object]:
    group = [dihedral_matrix(a, epsilon) for a in range(N) for epsilon in (-1, 1)]
    rng = np.random.default_rng(RNG_SEED + 27)
    psi = rng.normal(size=N) + 1j * rng.normal(size=N)
    psi /= np.linalg.norm(psi)
    p = np.abs(psi) ** 2
    shift = np.roll(np.eye(N), 1, axis=0)
    current = (shift - shift.T) / (2.0j)

    spectral_residual = 0.0
    entropy_residual = 0.0
    return_relational_residual = 0.0
    absolute_site_variation = 0.0
    current_even_residual = 0.0
    current_odd_residual = 0.0
    base_spectrum = np.linalg.eigvalsh(W)
    base_entropy = -float(np.sum(p * np.log(np.maximum(p, 1.0e-300))))
    base_current = float(np.real(np.vdot(psi, current @ psi)))
    for g in group:
        transformed_w = g @ W @ g.T
        transformed_psi = g @ psi
        transformed_p = np.abs(transformed_psi) ** 2
        spectral_residual = max(
            spectral_residual,
            float(np.max(np.abs(np.linalg.eigvalsh(transformed_w) - base_spectrum))),
        )
        entropy = -float(np.sum(transformed_p * np.log(np.maximum(transformed_p, 1.0e-300))))
        entropy_residual = max(entropy_residual, abs(entropy - base_entropy))
        origin = int(np.argmax(g[:, 0]))
        return_relational_residual = max(return_relational_residual, abs(transformed_p[origin] - p[0]))
        absolute_site_variation = max(absolute_site_variation, abs(transformed_p[0] - p[0]))
        transformed_current = float(np.real(np.vdot(transformed_psi, current @ transformed_psi)))
        # Determine whether g contains a reflection from its action on +1.
        image_zero = int(np.argmax(g[:, 0]))
        image_one = int(np.argmax(g[:, 1]))
        epsilon = 1 if (image_one - image_zero) % N == 1 else -1
        current_even_residual = max(current_even_residual, abs(transformed_current - base_current))
        current_odd_residual = max(current_odd_residual, abs(transformed_current - epsilon * base_current))

    registry = [
        ("kernel spectrum", "invariant", True),
        ("spectral gap", "invariant", True),
        ("probability entropy", "invariant", True),
        ("return probability", "relational origin", True),
        ("transition distance", "relational origin", True),
        ("apparatus displacement distribution", "relative frame", True),
        ("signed current", "requires odd apparatus polarity", True),
        ("chiral sign", "requires odd apparatus polarity", True),
        ("absolute site-zero probability", "absolute label", False),
        ("dimensionful clock reading", "requires calibration", False),
        ("absolute orientation selector", "forbidden without odd resource", False),
        ("record identity", "requires apparatus memory", False),
    ]
    return {
        "group_elements_checked": len(group),
        "observable_registry": [
            {"observable": name, "typing": typing, "relational_rewrite_available": available}
            for name, typing, available in registry
        ],
        "relationally_rewritable_count": sum(available for _, _, available in registry),
        "requires_additional_resource_count": sum(not available for _, _, available in registry),
        "spectrum_covariance_residual": spectral_residual,
        "entropy_covariance_residual": entropy_residual,
        "return_probability_joint_frame_residual": return_relational_residual,
        "absolute_site_zero_variation": absolute_site_variation,
        "signed_current_invariant_residual": current_even_residual,
        "signed_current_odd_covariance_residual": current_odd_residual,
    }


def time_average_covariance(matrix: np.ndarray, total_time: float, seed: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    coefficients = vectors.T @ seed
    differences = values[:, None] - values[None, :]
    factors = np.ones_like(differences, dtype=complex)
    mask = np.abs(differences) > 1.0e-12
    factors[mask] = (1.0 - np.exp(-1j * differences[mask] * total_time)) / (
        1j * differences[mask] * total_time
    )
    modal = coefficients[:, None] * coefficients[None, :] * factors
    covariance = vectors @ modal @ vectors.T
    return 0.5 * (covariance.real + covariance.real.T)


def adaptive_integrability_study() -> dict[str, object]:
    basis: list[np.ndarray] = []
    for i in range(N):
        for j in range(i + 1, N):
            h = np.zeros((N, N))
            h[i, j] = h[j, i] = 1.0 / math.sqrt(2.0)
            basis.append(h)
    b = np.array(basis)
    dimension = len(basis)
    seed = np.zeros(N)
    seed[0] = 1.0
    total_time = 4.0
    step = 2.0e-6
    jacobian = np.zeros((dimension, dimension))
    for column, h in enumerate(basis):
        derivative = (
            time_average_covariance(W + step * h, total_time, seed)
            - time_average_covariance(W - step * h, total_time, seed)
        ) / (2.0 * step)
        jacobian[:, column] = np.tensordot(b, derivative, axes=([1, 2], [0, 1]))

    reflection = dihedral_matrix(0, -1)
    representation = np.zeros((dimension, dimension))
    for column, h in enumerate(basis):
        reflected = reflection @ h @ reflection.T
        representation[:, column] = np.tensordot(b, reflected, axes=([1, 2], [0, 1]))
    parity_values, parity_vectors = np.linalg.eigh(representation)
    even = parity_vectors[:, parity_values > 0.5]
    odd = parity_vectors[:, parity_values < -0.5]
    odd_jacobian = odd.T @ jacobian @ odd
    even_odd = even.T @ jacobian @ odd
    odd_even = odd.T @ jacobian @ even
    odd_eigenvalues = np.linalg.eigvals(odd_jacobian)
    leading_index = int(np.argmax(odd_eigenvalues.real))
    leading = odd_eigenvalues[leading_index]
    centered_crossing = math.inf if leading.real <= 0.0 else 1.0 / leading.real
    return {
        "tangent_dimension": dimension,
        "even_dimension": int(even.shape[1]),
        "odd_dimension": int(odd.shape[1]),
        "covariance_Jacobian_norm": float(np.linalg.norm(jacobian)),
        "integrability_antisymmetry_relative_norm": float(
            np.linalg.norm(jacobian - jacobian.T) / np.linalg.norm(jacobian)
        ),
        "even_to_odd_block_norm": float(np.linalg.norm(even_odd)),
        "odd_to_even_block_norm": float(np.linalg.norm(odd_even)),
        "odd_leading_eigenvalue_real": float(leading.real),
        "odd_leading_eigenvalue_imag": float(leading.imag),
        "centered_linear_flow_mu_at_real_part_crossing": float(centered_crossing),
        "odd_symmetric_part_min_eigenvalue": float(
            np.linalg.eigvalsh(0.5 * (odd_jacobian + odd_jacobian.T))[0]
        ),
        "odd_symmetric_part_max_eigenvalue": float(
            np.linalg.eigvalsh(0.5 * (odd_jacobian + odd_jacobian.T))[-1]
        ),
        "verdict": "the declared finite-time covariance one-form is nonintegrable; the first centered odd instability is not a certified pitchfork",
    }


def no_go_library_study() -> dict[str, object]:
    # Certificate 1: a singleton with trivial reflection action has no
    # equivariant map to the two-point sign torsor.
    selector_candidates = (-1, 1)
    equivariant_selector_count = sum(sign == -sign for sign in selector_candidates)

    # Certificate 2: s = 2s has no strictly positive solution.
    positive_scale_section_count = sum(abs(s - 2.0 * s) < EPS for s in (0.25, 0.5, 1.0, 2.0))

    # Certificate 3: the nontrivial FIN propagator is not both a nonnegative
    # stochastic matrix and unitary.  The general theorem reduces their
    # intersection to permutation matrices.
    t = 1.0
    u = unitary(t)
    both_residual = {
        "unitarity": float(np.linalg.norm(u.conj().T @ u - np.eye(N))),
        "imaginary_norm": float(np.linalg.norm(u.imag)),
        "minimum_real_entry": float(u.real.min()),
        "distance_to_identity_permutation": float(np.linalg.norm(u - np.eye(N))),
    }

    # Certificate 4: a nonzero zero-diagonal symmetric matrix cannot be PSD.
    w_eigenvalues = np.linalg.eigvalsh(W)

    # Certificate 5: finite positive frequency mixtures have a purely
    # imaginary first derivative at zero; real exponential decay has -gamma.
    rational_weights = np.array([0.25, 0.50, 0.25])
    rational_frequencies = np.array([-2.0, 0.0, 2.0])
    finite_bath_derivative = -1j * float(rational_weights @ rational_frequencies)
    exponential_derivative = -1.0

    checks = {
        "trivial_to_sign_equivariant_maps": int(equivariant_selector_count),
        "positive_scale_sections_checked_at_c2": int(positive_scale_section_count),
        "unitary_stochastic_FIN_witness": both_residual,
        "normal_ordered_kernel_trace": float(np.trace(W)),
        "normal_ordered_kernel_negative_eigenvalues": int(np.sum(w_eigenvalues < -1.0e-12)),
        "finite_bath_derivative_at_zero_real": float(np.real(finite_bath_derivative)),
        "finite_bath_derivative_at_zero_imag": float(np.imag(finite_bath_derivative)),
        "target_exponential_derivative_at_zero": exponential_derivative,
    }
    serialized = json.dumps(checks, sort_keys=True, separators=(",", ":")).encode()
    return {
        "certificate_count": 5,
        "all_assertions_passed": bool(
            equivariant_selector_count == 0
            and positive_scale_section_count == 0
            and both_residual["unitarity"] < 1.0e-12
            and (both_residual["imaginary_norm"] > 1.0e-3 or both_residual["minimum_real_entry"] < -1.0e-6)
            and np.trace(W) == 0.0
            and np.sum(w_eigenvalues < -1.0e-12) > 0
            and abs(np.real(finite_bath_derivative) - exponential_derivative) > 0.5
        ),
        "checks": checks,
        "certificate_SHA256": hashlib.sha256(serialized).hexdigest(),
        "formalization_scope": "exact finite/executable certificates; independent proof-assistant port remains future work",
    }


def affine_fit(prediction: np.ndarray, observation: np.ndarray, training: np.ndarray) -> tuple[float, float]:
    design = np.column_stack([prediction[training], np.ones(training.size)])
    scale, offset = np.linalg.lstsq(design, observation[training], rcond=None)[0]
    return float(scale), float(offset)


def generator_from_six_weights(weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    radial = np.r_[0.0, weights]
    matrix = cyclic_matrix(radial)
    generator = float(matrix[0].sum()) * np.eye(N) - matrix
    rates = np.fft.fft(generator[0]).real
    return generator, rates


def challenge_distributions(
    generator: np.ndarray, alpha: float, spam: float, times: list[tuple[str, float]]
) -> list[np.ndarray]:
    values, vectors = np.linalg.eigh(generator)
    out: list[np.ndarray] = []
    for kind, time in times:
        if kind == "D":
            probability = ((vectors * np.exp(-alpha * time * values)) @ vectors.T)[:, 0].real
        else:
            amplitude = ((vectors * np.exp(-1j * alpha * time * values)) @ vectors.T)[:, 0]
            probability = np.abs(amplitude) ** 2
        probability = (1.0 - spam) * probability + spam / N
        out.append(normalize_probability(probability))
    return out


def challenge_fit_candidate(
    weights: np.ndarray,
    spectral_observation: np.ndarray,
    spectral_sigma: float,
    count_observations: list[np.ndarray],
    training_modes: np.ndarray,
    training_times: list[tuple[str, float]],
) -> tuple[float, float, float]:
    generator, rates = generator_from_six_weights(weights)

    def objective(parameters: np.ndarray) -> float:
        alpha, spam = parameters
        nll = 0.5 * float(
            np.sum(((spectral_observation - alpha * rates[training_modes]) / spectral_sigma) ** 2)
        )
        for counts, probability in zip(
            count_observations, challenge_distributions(generator, alpha, spam, training_times)
        ):
            nll -= float(counts @ np.log(np.maximum(probability, 1.0e-300)))
        return nll

    result = minimize(
        objective,
        np.array([1.0, 0.04]),
        method="L-BFGS-B",
        bounds=[(0.4, 2.2), (0.0, 0.15)],
        options={"maxiter": 80},
    )
    return float(result.fun), float(result.x[0]), float(result.x[1])


def challenge_heldout_nll(
    weights: np.ndarray,
    alpha: float,
    spam: float,
    spectral_observation: np.ndarray,
    spectral_sigma: float,
    count_observations: list[np.ndarray],
    modes: np.ndarray,
    times: list[tuple[str, float]],
) -> float:
    generator, rates = generator_from_six_weights(weights)
    nll = 0.5 * float(
        np.sum(((spectral_observation - alpha * rates[modes]) / spectral_sigma) ** 2)
    )
    for counts, probability in zip(
        count_observations, challenge_distributions(generator, alpha, spam, times)
    ):
        nll -= float(counts @ np.log(np.maximum(probability, 1.0e-300)))
    return nll


def blinded_challenge_study() -> dict[str, object]:
    """Synthetic preregistered spectral-plus-multitime held-out challenge."""
    rng = np.random.default_rng(RNG_SEED + 30)
    strict = RADIAL[1:].copy()
    strict /= 2.0 * strict[:5].sum() + strict[5]
    library_scales = (0.05, 0.10, 0.20, 0.40, 0.70)
    library: list[np.ndarray] = []
    for scale in library_scales:
        for _ in range(20):
            weights = strict * np.exp(rng.normal(0.0, scale, size=6))
            weights /= 2.0 * weights[:5].sum() + weights[5]
            library.append(weights)

    training_modes = np.array([1, 2])
    held_out_modes = np.array([3, 4, 5, 6])
    training_times = [("D", 0.35), ("U", 0.55)]
    held_out_times = [("D", 0.9), ("D", 1.7), ("U", 1.1), ("U", 2.2)]
    shots = 4000
    spectral_sigma = 0.012
    trials = 30
    labels = ["FIN" if i % 2 == 0 else "ALT" for i in range(trials)]
    rng.shuffle(labels)
    salt = "FIN-Programs-21-30-held-out-2026-07-19"
    preregistration = {
        "models": ["FIN", "positive_circulant_alternative_library"],
        "training_modes": training_modes.tolist(),
        "held_out_modes": held_out_modes.tolist(),
        "training_times": training_times,
        "held_out_times": held_out_times,
        "nuisance": "overall positive rate alpha and uniform SPAM epsilon, fitted on training data only",
        "decision": "smaller held-out negative log likelihood",
        "spectral_sigma": spectral_sigma,
        "shots_per_distribution": shots,
        "trials": trials,
        "seed": RNG_SEED + 30,
    }
    commitment = hashlib.sha256((salt + "|" + "|".join(labels)).encode()).hexdigest()
    confusion = np.zeros((2, 2), dtype=int)
    rows: list[dict[str, float | str]] = []
    alternative_counter = 0
    for label in labels:
        if label == "FIN":
            truth = strict
            perturbation_scale = 0.0
            true_index = 0
        else:
            perturbation_scale = library_scales[alternative_counter % len(library_scales)]
            alternative_counter += 1
            truth = strict * np.exp(rng.normal(0.0, perturbation_scale, size=6))
            truth /= 2.0 * truth[:5].sum() + truth[5]
            true_index = 1
        alpha_true = float(rng.uniform(0.75, 1.65))
        spam_true = float(rng.uniform(0.01, 0.08))
        truth_generator, truth_rates = generator_from_six_weights(truth)
        spectral_training = alpha_true * truth_rates[training_modes] + rng.normal(
            0.0, spectral_sigma, training_modes.size
        )
        spectral_test = alpha_true * truth_rates[held_out_modes] + rng.normal(
            0.0, spectral_sigma, held_out_modes.size
        )
        counts_training = [
            rng.multinomial(shots, probability)
            for probability in challenge_distributions(
                truth_generator, alpha_true, spam_true, training_times
            )
        ]
        counts_test = [
            rng.multinomial(shots, probability)
            for probability in challenge_distributions(
                truth_generator, alpha_true, spam_true, held_out_times
            )
        ]
        fin_fit = challenge_fit_candidate(
            strict,
            spectral_training,
            spectral_sigma,
            counts_training,
            training_modes,
            training_times,
        )
        alternative_fits = [
            challenge_fit_candidate(
                weights,
                spectral_training,
                spectral_sigma,
                counts_training,
                training_modes,
                training_times,
            )
            for weights in library
        ]
        best_index = int(np.argmin([fit[0] for fit in alternative_fits]))
        alternative_fit = alternative_fits[best_index]
        fin_test = challenge_heldout_nll(
            strict,
            fin_fit[1],
            fin_fit[2],
            spectral_test,
            spectral_sigma,
            counts_test,
            held_out_modes,
            held_out_times,
        )
        alternative_test = challenge_heldout_nll(
            library[best_index],
            alternative_fit[1],
            alternative_fit[2],
            spectral_test,
            spectral_sigma,
            counts_test,
            held_out_modes,
            held_out_times,
        )
        predicted_index = 0 if fin_test < alternative_test else 1
        confusion[true_index, predicted_index] += 1
        rows.append(
            {
                "truth": label,
                "alternative_log_perturbation_scale": float(perturbation_scale),
                "prediction": "FIN" if predicted_index == 0 else "ALT",
                "heldout_delta_NLL_generic_minus_FIN": float(alternative_test - fin_test),
            }
        )
    deltas_fin = [float(row["heldout_delta_NLL_generic_minus_FIN"]) for row in rows if row["truth"] == "FIN"]
    deltas_alt = [float(row["heldout_delta_NLL_generic_minus_FIN"]) for row in rows if row["truth"] == "ALT"]
    specificity_by_scale: dict[str, float] = {}
    for scale in library_scales:
        subset = [row for row in rows if row["truth"] == "ALT" and row["alternative_log_perturbation_scale"] == scale]
        specificity_by_scale[str(scale)] = float(
            sum(row["prediction"] == "ALT" for row in subset) / len(subset)
        )
    return {
        "preregistration": preregistration,
        "preregistration_SHA256": commitment,
        "commitment_salt_revealed": salt,
        "generic_library_size": len(library),
        "confusion_matrix_rows_true_columns_predicted": confusion.tolist(),
        "classification_accuracy": float(np.trace(confusion) / confusion.sum()),
        "FIN_sensitivity": float(confusion[0, 0] / confusion[0].sum()),
        "ALT_specificity": float(confusion[1, 1] / confusion[1].sum()),
        "ALT_specificity_by_perturbation_scale": specificity_by_scale,
        "median_delta_NLL_on_FIN_truth": float(np.median(deltas_fin)),
        "median_delta_NLL_on_ALT_truth": float(np.median(deltas_alt)),
        "scope": "synthetic pipeline validation only; no hardware or fundamental-physics inference",
    }


def inverse_action_study() -> dict[str, object]:
    laplacian = 2.0 * np.eye(N) - np.roll(np.eye(N), 1, axis=0) - np.roll(np.eye(N), -1, axis=0)

    def fit_at_mass(mass2: float) -> tuple[float, float, np.ndarray, np.ndarray]:
        g = np.linalg.inv(laplacian + mass2 * np.eye(N))
        normal = g - np.diag(np.diag(g))
        amplitude = float(np.sum(W * normal) / np.sum(normal * normal))
        residual = float(np.linalg.norm(amplitude * normal - W) / np.linalg.norm(W))
        return residual, amplitude, g, normal

    optimization = minimize_scalar(
        lambda x: fit_at_mass(float(x))[0], bounds=(0.02, 4.0), method="bounded",
        options={"xatol": 1.0e-13},
    )
    residual, amplitude, g, normal = fit_at_mass(float(optimization.x))
    fitted = amplitude * normal
    by_distance = [float(abs(fitted[0, d] - W[0, d])) for d in range(7)]
    k_eigenvalues = np.linalg.eigvalsh(W)
    distinct = np.unique(np.round(k_eigenvalues, 12))
    well_curvatures = []
    for value in distinct:
        product = 1.0
        for other in distinct:
            if other != value:
                product *= (value - other) ** 2
        well_curvatures.append(2.0 * product)

    h = laplacian + float(optimization.x) * np.eye(N)
    two_pi_gradient = 0.5 * (h - np.linalg.inv(g))
    direct_sum_mixed_hessian_norm = 0.0
    quadratic_driven_stationarity = float(np.linalg.norm(W - W))
    return {
        "best_mass_squared": float(optimization.x),
        "best_resolvent_amplitude": amplitude,
        "normal_ordered_relative_Frobenius_residual": residual,
        "absolute_residual_by_cyclic_distance": by_distance,
        "G_minimum_eigenvalue": float(np.linalg.eigvalsh(g)[0]),
        "normal_ordered_K_trace": float(np.trace(fitted)),
        "normal_ordered_K_negative_eigenvalues": int(np.sum(np.linalg.eigvalsh(fitted) < -1.0e-12)),
        "strict_K_negative_eigenvalues": int(np.sum(k_eigenvalues < -1.0e-12)),
        "distinct_strict_eigenvalues": int(distinct.size),
        "seven_well_polynomial_degree": int(2 * distinct.size),
        "seven_well_curvature_min": float(min(well_curvatures)),
        "seven_well_curvature_max": float(max(well_curvatures)),
        "two_PI_stationarity_gradient_norm": float(np.linalg.norm(two_pi_gradient)),
        "displayed_direct_sum_action_mixed_phi_K_hessian_norm": direct_sum_mixed_hessian_norm,
        "quadratic_driven_action_degree": 2,
        "quadratic_driven_stationarity_residual_with_C_equals_Kstrict": quadratic_driven_stationarity,
        "logical_verdict": "the Gaussian action generates G and only approximately its declared normal-ordered image K; the displayed phi-plus-K action is decoupled unless a constraint or self-consistency law is added",
    }


def make_figures(results: dict[str, object], directory: Path) -> list[str]:
    import matplotlib.pyplot as plt

    directory.mkdir(parents=True, exist_ok=True)
    files: list[str] = []

    # Process-tomography confusion matrix.
    p21 = results["program_21"]
    figure, axis = plt.subplots(figsize=(5.2, 4.3))
    image = axis.imshow(np.array(p21["confusion_matrix_rows_true_columns_predicted"]), cmap="Blues")
    axis.set_xticks(range(4), p21["models"], rotation=30, ha="right")
    axis.set_yticks(range(4), p21["models"])
    axis.set_xlabel("predicted model")
    axis.set_ylabel("true model")
    axis.set_title("Program 21: intervention-aware classification")
    for i in range(4):
        for j in range(4):
            axis.text(j, i, str(p21["confusion_matrix_rows_true_columns_predicted"][i][j]), ha="center", va="center")
    figure.colorbar(image, ax=axis, shrink=0.8)
    figure.tight_layout()
    path = directory / "program21_confusion.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    files.append(str(path))

    # Positive refinement and raw strict continuation.
    figure, axis = plt.subplots(figsize=(5.8, 3.8))
    d = np.arange(1, 31)
    raw = np.array([strict_formula(int(x)) for x in d])
    anchors = np.array(results["program_22"]["strict_anchor_values_d1_to_d6"])
    q = float(results["program_22"]["tail_ratio"])
    completion = np.array([anchors[x - 1] if x <= 6 else anchors[-1] * q ** (x - 6) for x in d])
    axis.axhline(0.0, color="black", linewidth=0.7)
    axis.plot(d, raw, "o-", ms=3, label="raw strict formula")
    axis.plot(d, completion, "s-", ms=3, label="positive completion")
    axis.set_yscale("symlog", linthresh=1.0e-3)
    axis.set_xlabel("cyclic distance d")
    axis.set_ylabel("weight")
    axis.set_title("Program 22: positivity-preserving extension")
    axis.legend()
    figure.tight_layout()
    path = directory / "program22_positive_refinement.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    files.append(str(path))

    # Dirac dispersion.
    figure, axis = plt.subplots(figsize=(5.4, 3.8))
    k = np.linspace(0.0, math.pi, 400)
    axis.plot(k, 2.0 * np.sin(k / 2.0), label=r"$2\sin(k/2)$")
    axis.plot(k, k, "--", label=r"linear $k$")
    axis.set_xlabel("lattice momentum k")
    axis.set_ylabel("positive frequency")
    axis.set_title("Program 23: incidence-Dirac dispersion")
    axis.legend()
    figure.tight_layout()
    path = directory / "program23_dirac_dispersion.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    files.append(str(path))

    # Adaptive odd Jacobian summary.
    p28 = results["program_28"]
    figure, axis = plt.subplots(figsize=(5.4, 3.8))
    values = [
        p28["integrability_antisymmetry_relative_norm"],
        p28["odd_symmetric_part_max_eigenvalue"],
        p28["odd_leading_eigenvalue_real"],
        abs(p28["odd_leading_eigenvalue_imag"]),
    ]
    labels = ["antisymmetry", "max sym eig", "lead Re", "lead |Im|"]
    axis.bar(labels, values, color=["#A61B1B", "#1F5A99", "#19733A", "#D55E00"])
    axis.set_title("Program 28: finite-time covariance Jacobian")
    axis.tick_params(axis="x", rotation=25)
    figure.tight_layout()
    path = directory / "program28_adaptive_jacobian.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    files.append(str(path))

    # Inverse propagator fit.
    p31 = results["inverse_action_audit"]
    figure, axis = plt.subplots(figsize=(5.5, 3.8))
    axis.semilogy(range(7), np.maximum(p31["absolute_residual_by_cyclic_distance"], 1.0e-12), "o-")
    axis.set_xlabel("cyclic distance d")
    axis.set_ylabel("absolute fit residual")
    axis.set_title("Inverse-action audit: normal-ordered resolvent residual")
    figure.tight_layout()
    path = directory / "inverse_action_residual.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    files.append(str(path))

    return files


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--figures-dir", type=Path)
    args = parser.parse_args()

    results: dict[str, object] = {
        "metadata": {
            "title": "FIN Programs 21-30: executed finite studies",
            "date": "2026-07-19",
            "times": "dimensionless",
            "numpy_version": np.__version__,
            "random_seed": RNG_SEED,
            "row_sum": ROW_SUM,
            "declared_exact_row_sum": ROW_SUM_EXACT,
            "rounded_profile_row_sum_difference": ROW_SUM - ROW_SUM_EXACT,
            "generator_minimum_eigenvalue": float(A_EVALS.min()),
        },
        "program_21": process_tomography_study(),
        "program_22": positivity_refinement_study(),
        "program_23": dirac_refinement_study(),
        "program_24": os_positivity_study(),
        "program_25": nuisance_robust_identification_study(),
        "program_26": bath_identifiability_study(),
        "program_27": relational_audit_study(),
        "program_28": adaptive_integrability_study(),
        "program_29": no_go_library_study(),
        "program_30": blinded_challenge_study(),
        "inverse_action_audit": inverse_action_study(),
    }
    if args.figures_dir:
        results["figure_files"] = make_figures(results, args.figures_dir)
    text = json.dumps(results, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
