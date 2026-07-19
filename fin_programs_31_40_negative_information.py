#!/usr/bin/env python3
"""Executed FIN Programs 31--40: negative information coupling.

This deterministic suite treats negative information coupling as a new,
explicitly declared operational hypothesis.  It does not assume a completed
legacy-to-strict bridge, transfer legacy physical roles, source a selector, or
derive physical units.  Every time and coupling in the calculations is
dimensionless.

Dependencies: NumPy, SciPy, Matplotlib.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import minimize, minimize_scalar


N = 12
ALPHA_GEO = 4.0 * math.log(2.0)
BETA_TORS = 0.01
OMEGA_LEGACY = math.pi / 4.0
PHI_LEGACY = math.pi / 6.0
OMEGA_STRICT = 0.18575
PHI_STRICT = 0.16250
BETA_STRICT = 1.0
ETA_STRICT = 1.8
SEED = 20260720


def legacy_weight(d: int | float) -> float:
    if float(d) == 0.0:
        return 0.0
    return ALPHA_GEO * math.cos(OMEGA_LEGACY * d + PHI_LEGACY) / (
        1.0 + BETA_TORS * d
    )


def strict_weight(d: int | float) -> float:
    if float(d) == 0.0:
        return 0.0
    return math.cos(OMEGA_STRICT * d + PHI_STRICT) / (
        1.0 + BETA_STRICT * d**ETA_STRICT
    )


def cyclic_distance(i: int, j: int, n: int = N) -> int:
    delta = abs(i - j)
    return min(delta, n - delta)


def cyclic_matrix_from_function(function, n: int = N) -> np.ndarray:
    return np.array(
        [
            [0.0 if i == j else function(cyclic_distance(i, j, n)) for j in range(n)]
            for i in range(n)
        ],
        dtype=float,
    )


def interval_matrix_from_function(function, n: int = N) -> np.ndarray:
    return np.array(
        [[0.0 if i == j else function(abs(i - j)) for j in range(n)] for i in range(n)],
        dtype=float,
    )


def laplacian(weights: np.ndarray) -> np.ndarray:
    return np.diag(weights.sum(axis=1)) - weights


W_LEGACY = cyclic_matrix_from_function(legacy_weight)
W_STRICT = cyclic_matrix_from_function(strict_weight)
A_LEGACY = laplacian(W_LEGACY)
A_STRICT = laplacian(W_STRICT)
UNIFORM = np.full(N, 1.0 / N)


def entropy(probability: np.ndarray) -> float:
    p = np.asarray(probability, dtype=float)
    p = p[p > 1.0e-300]
    return float(-np.sum(p * np.log(p)))


def relative_information(probability: np.ndarray) -> float:
    return math.log(N) - entropy(probability)


def trace_distance(rho: np.ndarray, sigma: np.ndarray) -> float:
    values = np.linalg.eigvalsh(0.5 * ((rho - sigma) + (rho - sigma).conj().T))
    return float(0.5 * np.sum(np.abs(values)))


def von_neumann_entropy(rho: np.ndarray) -> float:
    values = np.linalg.eigvalsh(0.5 * (rho + rho.conj().T)).real
    values = values[values > 1.0e-14]
    return float(-np.sum(values * np.log(values)))


def normalize_distribution(p: np.ndarray) -> np.ndarray:
    q = np.maximum(np.asarray(p, dtype=float).real, 0.0)
    return q / q.sum()


def program31_loss_only_bridge() -> dict[str, object]:
    distances = np.arange(1, 7, dtype=float)
    legacy = np.array([legacy_weight(d) for d in distances])
    strict = np.array([strict_weight(d) for d in distances])
    ratio = strict / legacy

    positive_mismatch = np.sign(legacy) != np.sign(strict)
    negative_mismatch = np.sign(-legacy) != np.sign(strict)
    # Independent nonnegative retentions can fit same-sign entries exactly and
    # must project opposite-sign targets to zero.
    positive_lower_bound = float(np.linalg.norm(strict[positive_mismatch]))
    negative_lower_bound = float(np.linalg.norm(strict[negative_mismatch]))
    strict_norm = float(np.linalg.norm(strict))

    # Loewner test on the canonical cycle and on an explicitly declared
    # interval extrapolation of the same analytic formulas.
    bridge_cycle = A_STRICT - A_LEGACY
    legacy_interval = laplacian(interval_matrix_from_function(legacy_weight))
    strict_interval = laplacian(interval_matrix_from_function(strict_weight))
    bridge_interval = strict_interval - legacy_interval
    eig_cycle = np.linalg.eigvalsh(bridge_cycle)
    eig_interval = np.linalg.eigvalsh(bridge_interval)

    return {
        "distances": distances.astype(int).tolist(),
        "legacy_row": legacy.tolist(),
        "strict_row": strict.tolist(),
        "strict_over_legacy": ratio.tolist(),
        "legacy_sign_pattern": "".join("+" if x > 0 else "-" for x in legacy),
        "strict_sign_pattern": "".join("+" if x > 0 else "-" for x in strict),
        "positive_normalization_sign_mismatches": int(np.sum(positive_mismatch)),
        "negative_normalization_sign_mismatches": int(np.sum(negative_mismatch)),
        "best_possible_positive_retention_residual_lower_bound": positive_lower_bound,
        "best_possible_positive_retention_relative_lower_bound": positive_lower_bound
        / strict_norm,
        "best_possible_negative_normalization_relative_lower_bound": negative_lower_bound
        / strict_norm,
        "cycle_bridge_minimum_eigenvalue": float(eig_cycle[0]),
        "cycle_bridge_maximum_eigenvalue": float(eig_cycle[-1]),
        "cycle_bridge_rank_tol_1e-10": int(np.sum(eig_cycle > 1.0e-10)),
        "interval_extrapolation_bridge_minimum_eigenvalue": float(eig_interval[0]),
        "interval_extrapolation_bridge_maximum_eigenvalue": float(eig_interval[-1]),
        "verdict": (
            "positive loss alone cannot map the legacy sign pattern to the strict row; "
            "the additive dissipative interpretation is positive on C12 but fails after "
            "the declared interval geometry extrapolation"
        ),
    }


def markov_diagnostics(weights: np.ndarray, times=(0.001, 0.01, 0.1, 1.0)) -> dict[str, object]:
    generator = laplacian(weights)
    eigenvalues = np.linalg.eigvalsh(generator)
    e0 = np.eye(N)[:, 0]
    time_rows = {}
    for t in times:
        transition = expm(-t * generator)
        p = transition @ e0
        time_rows[str(t)] = {
            "minimum_entry": float(np.min(transition)),
            "maximum_entry": float(np.max(transition)),
            "row_sum_residual": float(np.max(np.abs(transition.sum(axis=1) - 1.0))),
            "spectral_norm": float(np.linalg.norm(transition, 2)),
            "relative_information_if_probability": (
                relative_information(normalize_distribution(p))
                if float(np.min(transition)) >= -1.0e-12
                else None
            ),
        }
    positive = eigenvalues[eigenvalues > 1.0e-10]
    return {
        "row_sum": float(weights[0].sum()),
        "minimum_off_diagonal_weight": float(np.min(weights + np.eye(N) * 1.0e300)),
        "generator_eigenvalues": eigenvalues.tolist(),
        "minimum_generator_eigenvalue": float(eigenvalues[0]),
        "spectral_gap_if_psd": float(positive[0]) if eigenvalues[0] >= -1.0e-10 else None,
        "time_diagnostics": time_rows,
    }


def program32_legacy_markov() -> dict[str, object]:
    raw = markov_diagnostics(W_LEGACY)
    positive_weights = np.maximum(W_LEGACY, 0.0)
    absolute_weights = np.abs(W_LEGACY)
    positive = markov_diagnostics(positive_weights, times=(1.0,))
    absolute = markov_diagnostics(absolute_weights, times=(1.0,))

    e0 = np.eye(N)[:, 0]
    p_positive = expm(-laplacian(positive_weights)) @ e0
    p_absolute = expm(-laplacian(absolute_weights)) @ e0
    return {
        "raw_signed_legacy": raw,
        "positive_part_repair": {
            **positive,
            "relative_information_at_t1": relative_information(p_positive),
        },
        "absolute_value_repair": {
            **absolute,
            "relative_information_at_t1": relative_information(p_absolute),
        },
        "repair_gap_ratio_absolute_over_positive": float(
            absolute["spectral_gap_if_psd"] / positive["spectral_gap_if_psd"]
        ),
        "verdict": (
            "the raw signed legacy row is not a Markov generator; at least two "
            "inequivalent positivity repairs exist and predict very different mixing"
        ),
    }


def fit_retention_model(distances: np.ndarray, retention: np.ndarray, kind: str) -> dict[str, object]:
    if kind == "constant":
        objective = lambda gamma: np.linalg.norm(np.exp(-gamma * distances) - retention) / np.linalg.norm(retention)
        result = minimize_scalar(objective, bounds=(0.0, 5.0), method="bounded")
        prediction = np.exp(-result.x * distances)
        parameters = {"gamma": float(result.x)}
    elif kind == "stretched":
        def predict(x):
            gamma, power = np.exp(x)
            return np.exp(-gamma * distances**power)
        result = minimize(lambda x: np.linalg.norm(predict(x) - retention) / np.linalg.norm(retention), np.log([1.0, 0.7]), method="L-BFGS-B", bounds=[(-8, 4), (-4, 2)])
        prediction = predict(result.x)
        gamma, power = np.exp(result.x)
        parameters = {"gamma": float(gamma), "power": float(power)}
    elif kind == "power":
        def predict(x):
            scale, power = np.exp(x)
            return (1.0 + scale * distances) ** (-power)
        result = minimize(lambda x: np.linalg.norm(predict(x) - retention) / np.linalg.norm(retention), np.log([0.5, 2.0]), method="L-BFGS-B", bounds=[(-8, 4), (-4, 4)])
        prediction = predict(result.x)
        scale, power = np.exp(result.x)
        parameters = {"scale": float(scale), "power": float(power)}
    elif kind == "rational":
        def predict(x):
            beta, eta = np.exp(x)
            return 1.0 / (1.0 + beta * distances**eta)
        result = minimize(lambda x: np.linalg.norm(predict(x) - retention) / np.linalg.norm(retention), np.log([1.0, 1.8]), method="L-BFGS-B", bounds=[(-8, 4), (-4, 4)])
        prediction = predict(result.x)
        beta, eta = np.exp(result.x)
        parameters = {"beta": float(beta), "eta": float(eta)}
    else:
        raise ValueError(kind)
    return {
        "parameters": parameters,
        "relative_L2_error": float(np.linalg.norm(prediction - retention) / np.linalg.norm(retention)),
        "maximum_relative_error": float(np.max(np.abs(prediction - retention) / retention)),
        "prediction": prediction.tolist(),
    }


def program33_hazard_reconstruction() -> dict[str, object]:
    distances = np.arange(1, 13, dtype=float)
    legacy_envelope = 1.0 / (1.0 + BETA_TORS * distances)
    strict_envelope = 1.0 / (1.0 + distances**ETA_STRICT)
    retention = strict_envelope / legacy_envelope
    previous = np.concatenate([[1.0], retention[:-1]])
    local_retention = retention / previous
    hazard = -np.log(local_retention)

    fits = {
        kind: fit_retention_model(distances, retention, kind)
        for kind in ("constant", "stretched", "power", "rational")
    }

    rng = np.random.default_rng(SEED + 33)
    bootstrap_eta = []
    bootstrap_beta = []
    for _ in range(200):
        noisy = retention * np.exp(rng.normal(0.0, 0.01, size=retention.size))
        fit = fit_retention_model(distances, noisy, "rational")
        bootstrap_eta.append(fit["parameters"]["eta"])
        bootstrap_beta.append(fit["parameters"]["beta"])

    return {
        "distances": distances.astype(int).tolist(),
        "legacy_envelope": legacy_envelope.tolist(),
        "strict_envelope": strict_envelope.tolist(),
        "completion_retention": retention.tolist(),
        "local_retention": local_retention.tolist(),
        "discrete_hazard": hazard.tolist(),
        "all_hazards_nonnegative": bool(np.all(hazard >= -1.0e-14)),
        "candidate_fits": fits,
        "bootstrap_1pct_rational_eta_quantiles": np.quantile(bootstrap_eta, [0.025, 0.5, 0.975]).tolist(),
        "bootstrap_1pct_rational_beta_quantiles": np.quantile(bootstrap_beta, [0.025, 0.5, 0.975]).tolist(),
        "exact_declared_retention_formula": "(1+0.01 d)/(1+d^1.8)",
        "verdict": (
            "a positive cumulative-loss envelope exists exactly, but it is a "
            "reparameterization and does not select a microscopic loss mechanism or phase map"
        ),
    }


def feedback_weights(gamma: float, kind: str) -> np.ndarray:
    weights = np.zeros_like(W_STRICT)
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            d = cyclic_distance(i, j)
            base = strict_weight(d)
            if kind == "exponential":
                weights[i, j] = base * math.exp(-gamma * d)
            elif kind == "linear":
                weights[i, j] = base * (1.0 - gamma * d)
            else:
                raise ValueError(kind)
    return weights


def program34_stability_cone() -> dict[str, object]:
    selected = [0.0, 0.1, 1.0 / 6.0, 0.2, 0.4]
    rows: dict[str, dict[str, dict[str, object]]] = {"exponential": {}, "linear": {}}
    times = (0.001, 0.01, 0.1, 1.0)
    for kind in rows:
        for gamma in selected:
            weights = feedback_weights(gamma, kind)
            generator = laplacian(weights)
            eigenvalues = np.linalg.eigvalsh(generator)
            positives = eigenvalues[eigenvalues > 1.0e-10]
            transition_minima = {
                str(t): float(np.min(expm(-t * generator))) for t in times
            }
            rows[kind][f"{gamma:.12g}"] = {
                "minimum_off_diagonal_weight": float(np.min(weights + np.eye(N) * 1.0e300)),
                "minimum_generator_eigenvalue": float(eigenvalues[0]),
                "spectral_gap_first_positive": float(positives[0]) if positives.size else None,
                "transition_minimum_entries": transition_minima,
            }

    gamma_grid = np.linspace(0.0, 0.5, 101)
    phase_rows = []
    rng = np.random.default_rng(SEED + 34)
    kl_violations = {"exponential": 0, "linear_admissible": 0}
    for gamma in gamma_grid:
        record = {"gamma": float(gamma)}
        for kind in ("exponential", "linear"):
            weights = feedback_weights(gamma, kind)
            generator = laplacian(weights)
            eig = np.linalg.eigvalsh(generator)
            record[f"{kind}_minimum_weight"] = float(np.min(weights + np.eye(N) * 1.0e300))
            record[f"{kind}_minimum_eigenvalue"] = float(eig[0])
        phase_rows.append(record)

    for kind, gamma in (("exponential", 0.4), ("linear_admissible", 1.0 / 6.0)):
        actual_kind = "exponential" if kind == "exponential" else "linear"
        generator = laplacian(feedback_weights(gamma, actual_kind))
        for _ in range(100):
            p0 = rng.dirichlet(np.ones(N))
            infos = [relative_information(expm(-t * generator) @ p0) for t in np.linspace(0, 2, 21)]
            kl_violations[kind] += int(np.any(np.diff(infos) > 1.0e-10))

    return {
        "selected_rows": rows,
        "linear_exact_markov_threshold": 1.0 / 6.0,
        "phase_diagram": phase_rows,
        "relative_information_monotonicity_violations_over_100_trials": kl_violations,
        "verdict": (
            "positive exponential retention preserves the graph-Laplacian cone; "
            "linear subtraction ceases to be Markov-admissible above gamma=1/6"
        ),
    }


def program35_dilation_ledger() -> dict[str, object]:
    e0 = np.eye(N)[:, 0]
    times = np.linspace(0.0, 3.0, 61)
    curves = []
    for t in times:
        transition = expm(-t * A_STRICT)
        p = normalize_distribution(transition @ e0)
        h = entropy(p)
        curves.append(
            {
                "time": float(t),
                "system_entropy": h,
                "relative_information": math.log(N) - h,
                "system_purity": float(p @ p),
                "environment_entropy_pure_dilation": h,
                "system_environment_mutual_information": 2.0 * h,
            }
        )

    transition = expm(-A_STRICT)
    # Full isometry V|x> = sum_y sqrt(P(y|x)) |y>|x,y>.
    environment_dimension = N * N
    isometry = np.zeros((N * environment_dimension, N), dtype=complex)
    for x in range(N):
        for y in range(N):
            output = y * environment_dimension + x * N + y
            isometry[output, x] = math.sqrt(max(float(transition[y, x]), 0.0))
    isometry_residual = float(np.linalg.norm(isometry.conj().T @ isometry - np.eye(N)))
    psi = isometry @ e0
    psi_matrix = psi.reshape(N, environment_dimension)
    rho_system = psi_matrix @ psi_matrix.conj().T
    rho_environment = psi_matrix.conj().T @ psi_matrix
    p = normalize_distribution(transition @ e0)
    reduced_residual = float(np.linalg.norm(rho_system - np.diag(p)))
    h = entropy(p)

    return {
        "curves": curves,
        "t1": {
            "transition_minimum": float(np.min(transition)),
            "system_entropy": h,
            "relative_information": relative_information(p),
            "system_purity": float(p @ p),
            "joint_entropy": 0.0,
            "environment_entropy": von_neumann_entropy(rho_environment),
            "mutual_information": 2.0 * h,
        },
        "environment_dimension_of_declared_isometry": environment_dimension,
        "isometry_residual": isometry_residual,
        "reduced_channel_residual": reduced_residual,
        "verdict": (
            "the reduced system loses distinguishability while an explicit isometric "
            "dilation retains a pure joint state; fundamental deletion is not identifiable"
        ),
    }


SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
SIGMA_Y = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
SIGMA_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def partial_trace_environment(rho: np.ndarray, system_dimension: int = N) -> np.ndarray:
    tensor = rho.reshape(system_dimension, 2, system_dimension, 2)
    return np.trace(tensor, axis1=1, axis2=3)


def partial_trace_system(rho: np.ndarray, system_dimension: int = N) -> np.ndarray:
    tensor = rho.reshape(system_dimension, 2, system_dimension, 2)
    return np.trace(tensor, axis1=0, axis2=2)


def evolved_reduced(
    sign: float,
    time: float,
    system_state: np.ndarray,
    environment_state: np.ndarray,
    bridge: np.ndarray,
    coupling: float,
) -> tuple[np.ndarray, np.ndarray]:
    environment_hamiltonian = 0.7 * SIGMA_Z
    total = (
        np.kron(A_STRICT, I2)
        + np.kron(np.eye(N), environment_hamiltonian)
        + sign * coupling * np.kron(bridge, SIGMA_X)
    )
    unitary = expm(-1j * time * total)
    initial = np.kron(system_state, environment_state)
    final = unitary @ initial @ unitary.conj().T
    return partial_trace_environment(final), partial_trace_system(final)


def program36_sign_gauge() -> dict[str, object]:
    bridge = A_STRICT - A_LEGACY
    bridge = bridge / np.linalg.norm(bridge, "fro")
    coupling = 0.8
    tau_symmetric = 0.5 * (I2 + 0.35 * SIGMA_Z)
    tau_broken = 0.5 * (I2 + 0.35 * SIGMA_X)

    e0 = np.eye(N)[:, 0]
    e1 = np.eye(N)[:, 1]
    plus = np.ones(N) / math.sqrt(N)
    pair = (e0 + 1j * e1) / math.sqrt(2.0)
    probes = [np.outer(v, v.conj()) for v in (e0, plus, pair)]
    times = np.linspace(0.0, 4.0, 41)
    curves = []
    for t in times:
        system_sym = []
        system_broken = []
        env_odd = []
        for rho in probes:
            plus_s, plus_e = evolved_reduced(+1.0, t, rho, tau_symmetric, bridge, coupling)
            minus_s, minus_e = evolved_reduced(-1.0, t, rho, tau_symmetric, bridge, coupling)
            system_sym.append(trace_distance(plus_s, minus_s))
            env_odd.append(abs(float(np.trace((plus_e - minus_e) @ SIGMA_X).real)))

            plus_b, _ = evolved_reduced(+1.0, t, rho, tau_broken, bridge, coupling)
            minus_b, _ = evolved_reduced(-1.0, t, rho, tau_broken, bridge, coupling)
            system_broken.append(trace_distance(plus_b, minus_b))
        curves.append(
            {
                "time": float(t),
                "maximum_system_trace_distance_symmetric_bath": float(max(system_sym)),
                "maximum_system_trace_distance_polarized_bath": float(max(system_broken)),
                "maximum_environment_odd_record_difference": float(max(env_odd)),
            }
        )

    return {
        "coupling_magnitude": coupling,
        "bridge_Frobenius_norm": float(np.linalg.norm(bridge, "fro")),
        "curves": curves,
        "maximum_system_sign_difference_symmetric_bath": max(
            row["maximum_system_trace_distance_symmetric_bath"] for row in curves
        ),
        "maximum_system_sign_difference_polarized_bath": max(
            row["maximum_system_trace_distance_polarized_bath"] for row in curves
        ),
        "maximum_environment_odd_record_difference": max(
            row["maximum_environment_odd_record_difference"] for row in curves
        ),
        "verdict": (
            "the coupling sign is an exact system-only gauge for a parity-symmetric bath; "
            "a calibrated odd environment record or parity-breaking preparation is required"
        ),
    }


def controlled_rate_matrix(theta: float, field: np.ndarray) -> np.ndarray:
    rates = np.zeros((N, N), dtype=float)
    for destination in range(N):
        for source in range(N):
            if destination == source:
                continue
            base = W_STRICT[destination, source]
            rates[destination, source] = base * math.exp(
                -0.5 * theta * (field[destination] - field[source])
            )
    for source in range(N):
        rates[source, source] = -np.sum(rates[:, source])
    return rates


def entropy_production_rate(p: np.ndarray, rates: np.ndarray) -> float:
    sigma = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            forward = rates[i, j] * p[j]
            reverse = rates[j, i] * p[i]
            if forward > 1.0e-300 and reverse > 1.0e-300:
                sigma += (forward - reverse) * math.log(forward / reverse)
    return float(sigma)


def program37_feedback_ledger() -> dict[str, object]:
    field = np.cos(2.0 * math.pi * np.arange(N) / N)
    kappa = 0.8
    theta0 = 1.4
    p0 = np.eye(N)[:, 0]

    def rhs(_time, state):
        p = state[:N]
        theta = state[-1]
        rates = controlled_rate_matrix(theta, field)
        return np.concatenate([rates @ p, [-kappa * theta]])

    times = np.linspace(0.0, 6.0, 121)
    solution = solve_ivp(rhs, (times[0], times[-1]), np.concatenate([p0, [theta0]]), t_eval=times, rtol=1.0e-10, atol=1.0e-12)
    curves = []
    minimum_sigma = math.inf
    for index, t in enumerate(times):
        p = normalize_distribution(solution.y[:N, index])
        theta = float(solution.y[-1, index])
        rates = controlled_rate_matrix(theta, field)
        stationary = np.exp(-theta * field)
        stationary /= stationary.sum()
        kl = float(np.sum(p * np.log(np.maximum(p, 1.0e-300) / stationary)))
        sigma = entropy_production_rate(p, rates)
        minimum_sigma = min(minimum_sigma, sigma)
        curves.append(
            {
                "time": float(t),
                "theta": theta,
                "controller_quadratic_loss": 0.5 * theta * theta,
                "relative_entropy_to_instantaneous_stationary_state": kl,
                "markov_entropy_production_rate": sigma,
            }
        )

    omega = 1.7
    jacobian = np.array([[-kappa, omega], [-omega, -kappa]])
    eigenvalues = np.linalg.eigvals(jacobian)
    integrability_defect = float(np.linalg.norm(jacobian - jacobian.T) / np.linalg.norm(jacobian))
    angles = np.linspace(0.0, 2.0 * math.pi, 4001)
    points = np.column_stack([np.cos(angles), np.sin(angles)])
    forces = points @ jacobian.T
    increments = np.diff(points, axis=0)
    midpoint_forces = 0.5 * (forces[:-1] + forces[1:])
    circulation = float(np.sum(midpoint_forces * increments))

    return {
        "curves": curves,
        "minimum_markov_entropy_production_rate": minimum_sigma,
        "controller_initial_loss": 0.5 * theta0 * theta0,
        "controller_final_loss": curves[-1]["controller_quadratic_loss"],
        "stable_nonconservative_feedback": {
            "jacobian": jacobian.tolist(),
            "eigenvalues_real_imag": [[float(z.real), float(z.imag)] for z in eigenvalues],
            "integrability_defect": integrability_defect,
            "unit_circle_work": circulation,
        },
        "verdict": (
            "negative feedback can be stable and entropy-production compatible, but stable "
            "feedback need not be variational and its target/sign are additional inputs"
        ),
    }


def operator_from_signed_feedback(g: float) -> np.ndarray:
    weights = np.zeros_like(W_STRICT)
    for i in range(N):
        for j in range(N):
            if i != j:
                d = cyclic_distance(i, j)
                weights[i, j] = strict_weight(d) * math.exp(g * d / 6.0)
    return laplacian(weights)


def program38_wave_diffusion_observer() -> dict[str, object]:
    e0 = np.eye(N)[:, 0]
    e3 = np.eye(N)[:, 3]
    small_times = np.geomspace(1.0e-4, 0.08, 30)
    rows = {}
    for g in (-0.5, 0.0, 0.5):
        operator = operator_from_signed_feedback(g)
        diff_survival = []
        wave_survival = []
        for t in small_times:
            diff_survival.append(float(expm(-t * operator)[0, 0]))
            wave_survival.append(float(abs(expm(-1j * t * operator)[0, 0]) ** 2))
        diff_escape = np.maximum(1.0 - np.array(diff_survival), 1.0e-300)
        wave_escape = np.maximum(1.0 - np.array(wave_survival), 1.0e-300)
        diff_exponent = float(np.polyfit(np.log(small_times[:18]), np.log(diff_escape[:18]), 1)[0])
        wave_exponent = float(np.polyfit(np.log(small_times[:18]), np.log(wave_escape[:18]), 1)[0])

        ck = {}
        for t in (0.2, 0.5, 1.0):
            markov = expm(-t * operator)
            unitary_population = np.abs(expm(-1j * t * operator)) ** 2
            ck[str(t)] = {
                "diffusion_CK_residual": float(np.linalg.norm(expm(-2.0 * t * operator) - markov @ markov)),
                "unitary_population_CK_residual": float(
                    np.linalg.norm(np.abs(expm(-2j * t * operator)) ** 2 - unitary_population @ unitary_population)
                ),
            }

        phases = np.linspace(0.0, 2.0 * math.pi, 181)
        probabilities = []
        unitary = expm(-2.0j * operator)
        for phase in phases:
            psi = (e0 + np.exp(1j * phase) * e3) / math.sqrt(2.0)
            probabilities.append(np.abs(unitary @ psi) ** 2)
        probabilities = np.asarray(probabilities)
        max_p = probabilities.max(axis=0)
        min_p = probabilities.min(axis=0)
        visibility = (max_p - min_p) / np.maximum(max_p + min_p, 1.0e-300)
        detector = int(np.argmax(visibility))
        rows[str(g)] = {
            "diffusive_escape_exponent": diff_exponent,
            "unitary_escape_exponent": wave_exponent,
            "diffusive_linear_escape_coefficient": float(operator[0, 0]),
            "unitary_quadratic_escape_coefficient": float((operator @ operator)[0, 0] - operator[0, 0] ** 2),
            "CK_tests": ck,
            "maximum_double_path_visibility": float(visibility[detector]),
            "best_detector": detector,
            "small_times": small_times.tolist(),
            "diffusion_survival": diff_survival,
            "unitary_survival": wave_survival,
            "phase_scan": phases.tolist(),
            "best_detector_probability": probabilities[:, detector].tolist(),
        }

    return {
        "feedback_values": [-0.5, 0.0, 0.5],
        "rows": rows,
        "verdict": (
            "the same signed spatial parameter changes both branches but does not select "
            "the temporal category; linear versus quadratic escape and CK tests remain decisive"
        ),
    }


def poisson_means(parameters: np.ndarray, reference: bool, polarity: float = 1.0) -> np.ndarray:
    g, clock, background = parameters
    b = np.array([1.0, 0.7, -0.4, 1.2, 0.2, -0.9])
    c = polarity * np.array([0.2, 1.0, 0.8, -0.5, 1.3, 0.4])
    h = 0.65 if reference else 0.0
    exposure = np.array([35.0, 42.0, 31.0, 46.0, 38.0, 44.0])
    return 7.0 + background + exposure * (clock * (g * b + h * c)) ** 2


def numerical_fisher(parameters: np.ndarray, reference: bool) -> np.ndarray:
    means = poisson_means(parameters, reference)
    jacobian = np.zeros((means.size, parameters.size))
    for index in range(parameters.size):
        step = 1.0e-6 * max(1.0, abs(parameters[index]))
        plus = parameters.copy()
        minus = parameters.copy()
        plus[index] += step
        minus[index] -= step
        jacobian[:, index] = (poisson_means(plus, reference) - poisson_means(minus, reference)) / (2.0 * step)
    return jacobian.T @ (jacobian / means[:, None])


def poisson_nll(counts: np.ndarray, means: np.ndarray) -> float:
    return float(np.sum(means - counts * np.log(np.maximum(means, 1.0e-300))))


def profile_curve(counts: np.ndarray, reference: bool, polarity: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    grid = np.linspace(-0.9, 0.9, 181)
    values = []
    for g in grid:
        result = minimize(
            lambda x: poisson_nll(counts, poisson_means(np.array([g, math.exp(x[0]), math.exp(x[1])]), reference, polarity)),
            np.log([1.0, 1.0]),
            method="L-BFGS-B",
            bounds=[(math.log(0.6), math.log(1.4)), (math.log(0.01), math.log(8.0))],
        )
        values.append(float(result.fun))
    values = np.asarray(values)
    return grid, values - values.min()


def program39_identifiability() -> dict[str, object]:
    true = np.array([-0.42, 1.06, 1.3])
    fisher_blind = numerical_fisher(true, False)
    fisher_reference = numerical_fisher(true, True)
    evals_blind = np.linalg.eigvalsh(fisher_blind)
    evals_reference = np.linalg.eigvalsh(fisher_reference)

    rng = np.random.default_rng(SEED + 39)
    counts_blind = rng.poisson(poisson_means(true, False))
    counts_reference = rng.poisson(poisson_means(true, True))
    grid, profile_blind = profile_curve(counts_blind, False)
    _, profile_reference = profile_curve(counts_reference, True)
    _, profile_reflected = profile_curve(counts_reference, True, polarity=-1.0)

    negative_index = int(np.argmin(np.abs(grid - true[0])))
    positive_index = int(np.argmin(np.abs(grid + true[0])))
    return {
        "true_parameters_g_clock_background": true.tolist(),
        "blind_fisher_eigenvalues": evals_blind.tolist(),
        "reference_fisher_eigenvalues": evals_reference.tolist(),
        "blind_fisher_rank_tol_1e-8": int(np.sum(evals_blind > 1.0e-8)),
        "reference_fisher_rank_tol_1e-8": int(np.sum(evals_reference > 1.0e-8)),
        "profile_grid": grid.tolist(),
        "blind_delta_NLL": profile_blind.tolist(),
        "reference_delta_NLL": profile_reference.tolist(),
        "reflected_apparatus_delta_NLL": profile_reflected.tolist(),
        "blind_sign_delta_NLL_at_true_magnitude": float(profile_blind[positive_index] - profile_blind[negative_index]),
        "reference_sign_delta_NLL_at_true_magnitude": float(profile_reference[positive_index] - profile_reference[negative_index]),
        "best_g_with_reference": float(grid[np.argmin(profile_reference)]),
        "best_g_after_unrecoded_apparatus_reflection": float(grid[np.argmin(profile_reflected)]),
        "verdict": (
            "a g-scaling null direction and global sign symmetry persist without a calibrated "
            "reference; reference interference restores local rank and sign sensitivity"
        ),
    }


def apply_spam(probability: np.ndarray, epsilon: float) -> np.ndarray:
    p = np.asarray(probability, dtype=float)
    return (1.0 - epsilon) * p + epsilon * np.full(p.size, 1.0 / p.size)


def challenge_blocks(model: str, clock: float, epsilon: float) -> tuple[list[np.ndarray], list[np.ndarray]]:
    g0 = 0.28
    if model in ("negative", "static_mimic"):
        g = -g0
    elif model == "positive":
        g = g0
    elif model == "zero":
        g = 0.0
    else:
        raise ValueError(model)
    operator = operator_from_signed_feedback(g)
    e0 = np.eye(N)[:, 0]
    e3 = np.eye(N)[:, 3]
    pair = (e0 + e3) / math.sqrt(2.0)

    def diffusion(t):
        return apply_spam(normalize_distribution(expm(-clock * t * operator) @ e0), epsilon)

    def wave(t, state=e0):
        return apply_spam(np.abs(expm(-1j * clock * t * operator) @ state) ** 2, epsilon)

    sign_signal = {"negative": 1.0, "positive": -1.0, "zero": 0.0, "static_mimic": 0.0}[model]

    def environment_record(t):
        amplitude = 0.075 * math.tanh(clock * t)
        p = 0.5 + sign_signal * amplitude
        p = (1.0 - epsilon) * p + 0.5 * epsilon
        return np.array([p, 1.0 - p])

    training = [diffusion(0.35), wave(0.55), environment_record(0.6)]
    held_out = [
        diffusion(1.05),
        wave(1.35),
        wave(1.8, pair),
        environment_record(1.4),
        environment_record(2.2),
    ]
    return training, held_out


def block_log_likelihood(counts: list[np.ndarray], probabilities: list[np.ndarray]) -> float:
    return float(
        sum(count @ np.log(np.maximum(probability, 1.0e-300)) for count, probability in zip(counts, probabilities))
    )


def program40_blinded_challenge() -> dict[str, object]:
    models = ("negative", "zero", "positive", "static_mimic")
    protocol = {
        "models": list(models),
        "g0": 0.28,
        "training": ["diffusion_t0.35", "unitary_t0.55", "environment_t0.6"],
        "held_out": ["diffusion_t1.05", "unitary_t1.35", "double_path_t1.8", "environment_t1.4", "environment_t2.2"],
        "clock_grid": [0.95, 0.975, 1.0, 1.025, 1.05],
        "spam_grid": [0.0, 0.01, 0.02, 0.03, 0.04, 0.05],
        "shots_per_block": 500,
        "trials_per_model": 60,
        "seed": SEED + 40,
        "decision": "maximum held-out likelihood after training-only nuisance fit",
        "scope": "synthetic pipeline test only",
    }
    protocol_text = json.dumps(protocol, sort_keys=True, separators=(",", ":"))
    protocol_hash = hashlib.sha256(protocol_text.encode("utf-8")).hexdigest()

    clock_grid = protocol["clock_grid"]
    spam_grid = protocol["spam_grid"]
    library = {
        model: {
            (clock, spam): challenge_blocks(model, clock, spam)
            for clock in clock_grid
            for spam in spam_grid
        }
        for model in models
    }
    rng = np.random.default_rng(protocol["seed"])
    confusion = np.zeros((len(models), len(models)), dtype=int)
    fitted_rows = []
    shots = protocol["shots_per_block"]
    for true_index, true_model in enumerate(models):
        for _ in range(protocol["trials_per_model"]):
            true_clock = float(rng.uniform(0.95, 1.05))
            true_spam = float(rng.uniform(0.0, 0.05))
            training_p, held_p = challenge_blocks(true_model, true_clock, true_spam)
            training_counts = [rng.multinomial(shots, p) for p in training_p]
            held_counts = [rng.multinomial(shots, p) for p in held_p]
            held_scores = []
            choices = []
            for candidate in models:
                best = None
                for nuisance, (candidate_training, candidate_held) in library[candidate].items():
                    score = block_log_likelihood(training_counts, candidate_training)
                    if best is None or score > best[0]:
                        best = (score, nuisance, candidate_held)
                assert best is not None
                held_scores.append(block_log_likelihood(held_counts, best[2]))
                choices.append(best[1])
            prediction = int(np.argmax(held_scores))
            confusion[true_index, prediction] += 1
            fitted_rows.append(
                {
                    "true_model": true_model,
                    "predicted_model": models[prediction],
                    "true_clock": true_clock,
                    "true_spam": true_spam,
                    "held_score_margin": float(np.sort(held_scores)[-1] - np.sort(held_scores)[-2]),
                    "winning_training_nuisance": list(choices[prediction]),
                }
            )

    negative_system_train, negative_system_held = challenge_blocks("negative", 1.0, 0.0)
    mimic_system_train, mimic_system_held = challenge_blocks("static_mimic", 1.0, 0.0)
    # Remove the environment blocks: indices 2 in training and 3,4 in held-out.
    system_residual = max(
        np.linalg.norm(a - b)
        for a, b in zip(
            negative_system_train[:2] + negative_system_held[:3],
            mimic_system_train[:2] + mimic_system_held[:3],
        )
    )

    per_model_accuracy = {
        model: float(confusion[i, i] / confusion[i].sum()) for i, model in enumerate(models)
    }
    return {
        "preregistration": protocol,
        "preregistration_SHA256": protocol_hash,
        "confusion_matrix_rows_true_columns_predicted": confusion.tolist(),
        "classification_accuracy": float(np.trace(confusion) / confusion.sum()),
        "per_model_accuracy": per_model_accuracy,
        "negative_vs_static_mimic_system_only_signature_residual": float(system_residual),
        "median_winning_held_score_margin": float(np.median([row["held_score_margin"] for row in fitted_rows])),
        "fitted_trial_rows": fitted_rows,
        "verdict": (
            "the synthetic operational record can test a calibrated sign, while system-only "
            "kernel dynamics cannot distinguish negative coupling from an identical static mimic"
        ),
    }


def historical_arithmetic_audit() -> dict[str, object]:
    return {
        "exp_minus_2p9_times_7": math.exp(-2.9 * 7.0),
        "exp_minus_2p9_times_12": math.exp(-2.9 * 12.0),
        "diagram_reported_values": {"d7": 9.0e-4, "d12": 6.0e-6},
        "exact_cosine_zero_formula": "d=4/3+4n",
        "integer_claim_2_5_8_11_is_exact": False,
        "path_count_exponent_with_alpha_0p6": 1.6 - 0.6,
        "path_amplitude_exponent_needed_for_inverse_d": 2.6,
        "operational_loss_fields_present_in_diagram": 0,
        "operational_loss_fields_tested": [
            "state_space",
            "time_law",
            "channel_or_semigroup",
            "environment",
            "instrument",
            "measurement_record",
            "information_functional",
        ],
        "verdict": "historical motivation, not a theorem-level operational information-loss derivation",
    }


def generate_figures(results: dict[str, object], directory: Path) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    directory.mkdir(parents=True, exist_ok=True)
    files = []

    p31 = results["program31"]
    d = np.array(p31["distances"])
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.axhline(0, color="black", linewidth=0.8)
    ax.plot(d, p31["legacy_row"], "o-", label="canonical legacy")
    ax.plot(d, p31["strict_row"], "s-", label="strict gate")
    ax.set(xlabel="cyclic distance d", ylabel="kernel weight", title="Program 31: loss-only bridge sign obstruction")
    ax.legend()
    ax.grid(alpha=0.25)
    path = directory / "program31_loss_bridge.png"
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(str(path))

    raw_eig = np.array(results["program32"]["raw_signed_legacy"]["generator_eigenvalues"])
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.axhline(0, color="black", linewidth=0.8)
    ax.plot(np.arange(N), raw_eig, "o-", label="raw signed legacy")
    ax.set(xlabel="ordered spectral index", ylabel="generator eigenvalue", title="Program 32: raw legacy generator is indefinite")
    ax.grid(alpha=0.25); ax.legend()
    path = directory / "program32_legacy_spectrum.png"
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(str(path))

    p33 = results["program33"]
    d = np.array(p33["distances"])
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    axes[0].semilogy(d, p33["completion_retention"], "o-", label="exact retention")
    for name in ("stretched", "power", "rational"):
        axes[0].semilogy(d, p33["candidate_fits"][name]["prediction"], "--", label=name)
    axes[0].set(xlabel="d", ylabel="retention", title="Envelope completion"); axes[0].legend(fontsize=8); axes[0].grid(alpha=0.25)
    axes[1].plot(d, p33["discrete_hazard"], "o-", color="#A61B1B")
    axes[1].set(xlabel="d", ylabel="local hazard", title="Positive discrete loss hazard"); axes[1].grid(alpha=0.25)
    path = directory / "program33_hazard.png"
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(str(path))

    phase = results["program34"]["phase_diagram"]
    gamma = np.array([row["gamma"] for row in phase])
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    axes[0].plot(gamma, [row["exponential_minimum_weight"] for row in phase], label="exponential")
    axes[0].plot(gamma, [row["linear_minimum_weight"] for row in phase], label="linear")
    axes[0].axhline(0, color="black", linewidth=0.8); axes[0].axvline(1/6, color="#A61B1B", linestyle="--")
    axes[0].set(xlabel="feedback gamma", ylabel="minimum rate", title="Markov-rate boundary"); axes[0].legend(); axes[0].grid(alpha=0.25)
    axes[1].plot(gamma, [row["exponential_minimum_eigenvalue"] for row in phase], label="exponential")
    axes[1].plot(gamma, [row["linear_minimum_eigenvalue"] for row in phase], label="linear")
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set(xlabel="feedback gamma", ylabel="minimum eigenvalue", title="PSD boundary"); axes[1].legend(); axes[1].grid(alpha=0.25)
    path = directory / "program34_stability_cone.png"
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(str(path))

    ledger = results["program35"]["curves"]
    t = [row["time"] for row in ledger]
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.plot(t, [row["relative_information"] for row in ledger], label="system information D(p||u)")
    ax.plot(t, [row["environment_entropy_pure_dilation"] for row in ledger], label="environment entropy")
    ax.plot(t, [row["system_environment_mutual_information"] for row in ledger], label="mutual information")
    ax.set(xlabel="dimensionless time", ylabel="nats", title="Program 35: reduced loss is compatible with transfer")
    ax.legend(); ax.grid(alpha=0.25)
    path = directory / "program35_information_ledger.png"
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(str(path))

    p36 = results["program36"]["curves"]
    t = [row["time"] for row in p36]
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.semilogy(t, np.maximum([row["maximum_system_trace_distance_symmetric_bath"] for row in p36], 1e-16), label="system, symmetric bath")
    ax.semilogy(t, np.maximum([row["maximum_system_trace_distance_polarized_bath"] for row in p36], 1e-16), label="system, polarized bath")
    ax.semilogy(t, np.maximum([row["maximum_environment_odd_record_difference"] for row in p36], 1e-16), label="odd environment record")
    ax.set(xlabel="dimensionless time", ylabel="sign-discrimination diagnostic", title="Program 36: coupling-sign gauge")
    ax.legend(); ax.grid(alpha=0.25)
    path = directory / "program36_sign_gauge.png"
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(str(path))

    p38 = results["program38"]["rows"]["-0.5"]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    axes[0].loglog(p38["small_times"], 1-np.array(p38["diffusion_survival"]), label="diffusion escape")
    axes[0].loglog(p38["small_times"], 1-np.array(p38["unitary_survival"]), label="wave escape")
    axes[0].set(xlabel="t", ylabel="1-survival", title="Linear versus quadratic escape"); axes[0].legend(); axes[0].grid(alpha=0.25)
    axes[1].plot(p38["phase_scan"], p38["best_detector_probability"])
    axes[1].set(xlabel="relative slit phase", ylabel="detector probability", title="Unitary phase-sensitive record"); axes[1].grid(alpha=0.25)
    path = directory / "program38_wave_diffusion.png"
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(str(path))

    p39 = results["program39"]
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.plot(p39["profile_grid"], p39["blind_delta_NLL"], label="no reference")
    ax.plot(p39["profile_grid"], p39["reference_delta_NLL"], label="calibrated reference")
    ax.set(xlabel="signed coupling g", ylabel="profile delta NLL", title="Program 39: sign and scale identifiability")
    ax.set_ylim(0, min(60, max(p39["reference_delta_NLL"]))); ax.legend(); ax.grid(alpha=0.25)
    path = directory / "program39_identifiability.png"
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(str(path))

    confusion = np.array(results["program40"]["confusion_matrix_rows_true_columns_predicted"])
    labels = results["program40"]["preregistration"]["models"]
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    image = ax.imshow(confusion, cmap="Blues")
    for i in range(confusion.shape[0]):
        for j in range(confusion.shape[1]):
            ax.text(j, i, str(confusion[i, j]), ha="center", va="center")
    ax.set_xticks(range(len(labels)), labels, rotation=25, ha="right")
    ax.set_yticks(range(len(labels)), labels)
    ax.set(xlabel="predicted model", ylabel="true model", title="Program 40: blinded synthetic challenge")
    fig.colorbar(image, ax=ax)
    path = directory / "program40_confusion.png"
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(str(path))

    return files


def run_all() -> dict[str, object]:
    results = {
        "metadata": {
            "date": "2026-07-19",
            "random_seed": SEED,
            "numpy_version": np.__version__,
            "n": N,
            "scope": (
                "dimensionless conditional hypothesis tests; no physical units, strict selector, "
                "legacy role transfer, completed bridge, or fundamental information destruction"
            ),
        },
        "historical_arithmetic_audit": historical_arithmetic_audit(),
        "program31": program31_loss_only_bridge(),
        "program32": program32_legacy_markov(),
        "program33": program33_hazard_reconstruction(),
        "program34": program34_stability_cone(),
        "program35": program35_dilation_ledger(),
        "program36": program36_sign_gauge(),
        "program37": program37_feedback_ledger(),
        "program38": program38_wave_diffusion_observer(),
        "program39": program39_identifiability(),
        "program40": program40_blinded_challenge(),
    }
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="FIN_Programs_31_40_Negative_Information_Results.json")
    parser.add_argument("--figures", default="FIN_Programs_31_40_Figures")
    arguments = parser.parse_args()

    results = run_all()
    figure_files = generate_figures(results, Path(arguments.figures))
    results["figure_files"] = figure_files
    Path(arguments.output).write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"output": arguments.output, "figures": figure_files}, indent=2))


if __name__ == "__main__":
    main()
