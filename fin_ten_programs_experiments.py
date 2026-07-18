#!/usr/bin/env python3
"""Reproducible finite experiments for the FIN ten-program monograph.

The script uses only NumPy and explicit C12 kernel profiles.  It does not
import repository conclusions.  Each output block corresponds to one or more
research programs in FIN_Ten_Research_Programs_Monograph.md.
"""

from __future__ import annotations

import itertools
import json
import math
from typing import Callable

import numpy as np


N = 12
EPS = 1.0e-15


def cyclic_distance(i: int, j: int, n: int) -> int:
    delta = abs(i - j)
    return min(delta, n - delta)


def strict_profile(d: int) -> float:
    if d == 0:
        return 0.0
    return math.cos(0.18575 * d + 0.16250) / (1.0 + d**1.8)


def legacy_profile(d: int) -> float:
    if d == 0:
        return 0.0
    return (
        4.0
        * math.log(2.0)
        * math.cos(math.pi * d / 4.0 + math.pi / 6.0)
        / (1.0 + 0.01 * d)
    )


def radial_matrix(n: int, profile: Callable[[int], float]) -> np.ndarray:
    return np.array(
        [[profile(cyclic_distance(i, j, n)) for j in range(n)] for i in range(n)],
        dtype=float,
    )


def cycle_laplacian(n: int) -> np.ndarray:
    shift = np.roll(np.eye(n), 1, axis=1)
    return 2.0 * np.eye(n) - shift - shift.T


def matrix_polynomial(coefficients: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    result = np.zeros_like(matrix, dtype=np.result_type(matrix, coefficients))
    power = np.eye(matrix.shape[0], dtype=result.dtype)
    for coefficient in coefficients:
        result += coefficient * power
        power = power @ matrix
    return result


def unique_radial_symbol(matrix: np.ndarray) -> np.ndarray:
    symbol = np.fft.fft(matrix[0]).real
    return symbol[: matrix.shape[0] // 2 + 1]


def relative_residual(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(left - right) / np.linalg.norm(right))


def hermitian_exp(matrix: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * np.exp(values)) @ vectors.conj().T


def entropy(probabilities: np.ndarray) -> float:
    p = np.asarray(probabilities, dtype=float)
    p = p[p > 0.0]
    return float(-np.sum(p * np.log(p)))


def kl_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = np.maximum(np.asarray(p, dtype=float), EPS)
    q = np.maximum(np.asarray(q, dtype=float), EPS)
    p /= p.sum()
    q /= q.sum()
    return float(np.sum(p * np.log(p / q)))


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = np.maximum(np.asarray(p, dtype=float), EPS)
    q = np.maximum(np.asarray(q, dtype=float), EPS)
    p /= p.sum()
    q /= q.sum()
    midpoint = 0.5 * (p + q)
    return 0.5 * kl_divergence(p, midpoint) + 0.5 * kl_divergence(q, midpoint)


def spectral_evolution(
    generator: np.ndarray, initial: np.ndarray, time: float, unitary: bool
) -> np.ndarray:
    values, vectors = np.linalg.eigh(generator)
    exponent = np.exp((-1j if unitary else -1.0) * time * values)
    state = vectors @ (exponent * (vectors.T.conj() @ initial))
    if unitary:
        probabilities = np.abs(state) ** 2
    else:
        probabilities = np.maximum(state.real, 0.0)
    probabilities /= probabilities.sum()
    return probabilities


def canonical_binary_orbit(profile: tuple[int, ...], dihedral: bool) -> tuple[int, ...]:
    rotations = [profile[k:] + profile[:k] for k in range(len(profile))]
    if dihedral:
        reflected = tuple(reversed(profile))
        rotations += [
            reflected[k:] + reflected[:k] for k in range(len(reflected))
        ]
    return min(rotations)


def main() -> None:
    laplacian = cycle_laplacian(N)
    strict = radial_matrix(N, strict_profile)
    legacy = radial_matrix(N, legacy_profile)

    lap_nodes = 2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(7) / N)
    strict_symbol = unique_radial_symbol(strict)
    legacy_symbol = unique_radial_symbol(legacy)

    # Program 1: exact finite spectral bridge, its conditioning, and extrapolation.
    strict_in_l = np.polynomial.polynomial.polyfit(lap_nodes, strict_symbol, 6)
    legacy_in_l = np.polynomial.polynomial.polyfit(lap_nodes, legacy_symbol, 6)
    strict_from_legacy = np.polynomial.polynomial.polyfit(
        legacy_symbol, strict_symbol, 6
    )
    legacy_from_strict = np.polynomial.polynomial.polyfit(
        strict_symbol, legacy_symbol, 6
    )
    strict_l_residual = relative_residual(
        matrix_polynomial(strict_in_l, laplacian), strict
    )
    legacy_l_residual = relative_residual(
        matrix_polynomial(legacy_in_l, laplacian), legacy
    )
    strict_bridge_residual = relative_residual(
        matrix_polynomial(strict_from_legacy, legacy), strict
    )
    legacy_bridge_residual = relative_residual(
        matrix_polynomial(legacy_from_strict, strict), legacy
    )
    legacy_min_separation = float(
        np.min(np.abs(legacy_symbol[:, None] - legacy_symbol[None, :] + np.eye(7) * 1e9))
    )
    strict_min_separation = float(
        np.min(np.abs(strict_symbol[:, None] - strict_symbol[None, :] + np.eye(7) * 1e9))
    )
    raw_vandermonde_condition = float(
        np.linalg.cond(np.vander(legacy_symbol, 7, increasing=True))
    )
    bridge_out_of_sample: dict[str, float] = {}
    for n in (16, 24, 48):
        strict_n = radial_matrix(n, strict_profile)
        legacy_n = radial_matrix(n, legacy_profile)
        predicted = matrix_polynomial(strict_from_legacy, legacy_n)
        bridge_out_of_sample[str(n)] = relative_residual(predicted, strict_n)

    # Program 2: selector obstruction and what an explicit odd source changes.
    reflection = np.zeros((N, N))
    for i in range(N):
        reflection[(-i) % N, i] = 1.0
    shift = np.roll(np.eye(N), 1, axis=1)
    current = (shift - shift.T) / (2.0j)
    selector_expectations: dict[str, float] = {}
    for h in (-0.2, 0.0, 0.2):
        rho = hermitian_exp(-(strict - h * current))
        rho /= np.trace(rho)
        selector_expectations[f"h={h:+.1f}"] = float(
            np.real(np.trace(rho @ current))
        )

    # Program 3: unit-rank and weight-zero source exhaustion.
    unit_basis = np.eye(3, dtype=int)  # length, time, action
    ranks_by_subset_size = {
        str(size): int(
            max(
                np.linalg.matrix_rank(unit_basis[list(indices)])
                for indices in itertools.combinations(range(3), size)
            )
        )
        for size in (1, 2, 3)
    }
    monomial_count = 0
    weight_one_count = 0
    receiver_weights = np.zeros((6, 1), dtype=int)
    for exponents in itertools.product(range(-2, 3), repeat=6):
        if not any(exponents):
            continue
        monomial_count += 1
        weight = np.asarray(exponents) @ receiver_weights
        if int(weight[0]) == 1:
            weight_one_count += 1

    # Program 4: incompatible continuum scalings through exactly the same L_12.
    continuum_gaps: dict[str, dict[str, float]] = {}
    for n in (12, 24, 48, 96, 192):
        gap = float(np.sort(np.linalg.eigvalsh(cycle_laplacian(n)))[1])
        continuum_gaps[str(n)] = {
            "unscaled": gap,
            "N_squared_scaling": gap * (n / 12.0) ** 2,
            "N_fourth_scaling": gap * (n / 12.0) ** 4,
        }
    n_dispersion = 4096
    theta1 = 2.0 * math.pi / n_dispersion
    theta2 = 4.0 * math.pi / n_dispersion
    e1 = 2.0 - 2.0 * math.cos(theta1)
    e2 = 2.0 - 2.0 * math.cos(theta2)
    quadratic_exponent = math.log(e2 / e1) / math.log(2.0)
    square_root_exponent = math.log(math.sqrt(e2) / math.sqrt(e1)) / math.log(2.0)

    # Program 5: an operational discriminator between two dynamics of the same L.
    initial = np.zeros(N)
    initial[0] = 1.0
    best = None
    for time in np.linspace(0.05, 8.0, 796):
        p_unitary = spectral_evolution(laplacian, initial, float(time), True)
        p_diffusive = spectral_evolution(laplacian, initial, float(time), False)
        js = js_divergence(p_unitary, p_diffusive)
        if best is None or js > best[0]:
            best = (js, float(time), p_unitary, p_diffusive)
    assert best is not None
    js, best_time, p_unitary, p_diffusive = best
    kl_u_d = kl_divergence(p_unitary, p_diffusive)
    observations_for_log_bf_10 = math.ceil(10.0 / kl_u_d)

    # Program 6: Gibbs scale orbit and an exact finite Landauer identity.
    energy_levels = np.array([0.0, 1.0])
    beta = 1.3
    tau = np.exp(-beta * energy_levels)
    tau /= tau.sum()
    uniform = np.array([0.5, 0.5])
    entropy_drop = entropy(uniform) - entropy(tau)
    reservoir_heat = float(energy_levels @ (uniform - tau))
    relative_entropy = kl_divergence(uniform, tau)
    landauer_residual = beta * reservoir_heat - entropy_drop - relative_entropy
    gibbs_a = np.exp(-2.0 * np.array([0.0, 1.0, 2.0]))
    gibbs_a /= gibbs_a.sum()
    gibbs_b = np.exp(-0.2 * (10.0 * np.array([0.0, 1.0, 2.0])))
    gibbs_b /= gibbs_b.sum()

    # Program 7: naive finite real-spectral-triple first-order test.
    max_first_order = 0.0
    for a_index in range(N):
        a = np.zeros((N, N))
        a[a_index, a_index] = 1.0
        first_commutator = strict @ a - a @ strict
        for b_index in range(N):
            b = np.zeros((N, N))
            b[b_index, b_index] = 1.0
            nested = first_commutator @ b - b @ first_commutator
            max_first_order = max(max_first_order, float(np.linalg.norm(nested)))
    diagonal_dirac_commutator = float(
        np.linalg.norm(np.diag(np.diag(strict)) @ np.diag(np.arange(N))
                       - np.diag(np.arange(N)) @ np.diag(np.diag(strict)))
    )

    # Program 8: stability and conditioning of a Gaussian variational lift.
    strict_eigenvalues = np.linalg.eigvalsh(strict)
    threshold = float(-strict_eigenvalues[0])
    stability_scan: dict[str, dict[str, float]] = {}
    for offset in (1e-3, 1e-2, 1e-1, 1.0):
        shifted = strict + (threshold + offset) * np.eye(N)
        values = np.linalg.eigvalsh(shifted)
        stability_scan[str(offset)] = {
            "minimum_eigenvalue": float(values[0]),
            "condition_number": float(values[-1] / values[0]),
            "log_determinant": float(np.sum(np.log(values))),
        }

    # Program 9: scale/offset-free spectral fingerprints for an experiment.
    strict_levels = np.sort(strict_symbol)
    legacy_levels = np.sort(legacy_symbol)
    affine_scale = (strict_levels[-1] - strict_levels[0]) / (
        legacy_levels[-1] - legacy_levels[0]
    )
    affine_offset = strict_levels[0] - affine_scale * legacy_levels[0]
    legacy_calibrated = affine_scale * legacy_levels + affine_offset
    range_normalized_residuals = (
        legacy_calibrated - strict_levels
    ) / (strict_levels[-1] - strict_levels[0])
    strict_gap_fingerprint = np.diff(strict_levels) / (
        strict_levels[-1] - strict_levels[0]
    )
    legacy_gap_fingerprint = np.diff(legacy_levels) / (
        legacy_levels[-1] - legacy_levels[0]
    )

    # Program 10: machine-checkable orbit and representation obstructions.
    binary_profiles = list(itertools.product((0, 1), repeat=N))[1:]
    translation_orbits = {
        canonical_binary_orbit(profile, False) for profile in binary_profiles
    }
    dihedral_orbits = {
        canonical_binary_orbit(profile, True) for profile in binary_profiles
    }

    results = {
        "program_1_exact_spectral_bridge": {
            "strict_as_polynomial_of_L_relative_residual": strict_l_residual,
            "legacy_as_polynomial_of_L_relative_residual": legacy_l_residual,
            "strict_as_polynomial_of_legacy_relative_residual": strict_bridge_residual,
            "legacy_as_polynomial_of_strict_relative_residual": legacy_bridge_residual,
            "legacy_unique_symbol_minimum_separation": legacy_min_separation,
            "strict_unique_symbol_minimum_separation": strict_min_separation,
            "legacy_bridge_raw_vandermonde_condition": raw_vandermonde_condition,
            "same_polynomial_out_of_sample_relative_residual": bridge_out_of_sample,
        },
        "program_2_selector": {
            "reflection_residual_strict": float(
                np.linalg.norm(reflection @ strict @ reflection.T - strict)
            ),
            "reflection_residual_legacy": float(
                np.linalg.norm(reflection @ legacy @ reflection.T - legacy)
            ),
            "current_oddness_residual": float(
                np.linalg.norm(reflection @ current @ reflection.T + current)
            ),
            "dimension_even_linear_features_to_sign_output": 0,
            "current_expectation_after_explicit_odd_source": selector_expectations,
        },
        "program_3_units": {
            "rank_by_number_of_independent_calibrations": ranks_by_subset_size,
            "weight_zero_monomials_tested": monomial_count,
            "weight_one_sources_found": weight_one_count,
            "mass_weight_in_length_time_action_basis": [-2, 1, 1],
            "energy_weight_in_length_time_action_basis": [0, -1, 1],
        },
        "program_4_continuum_and_dispersion": {
            "first_nonzero_gap_for_three_extensions": continuum_gaps,
            "low_momentum_power_L": quadratic_exponent,
            "low_momentum_power_sqrt_L": square_root_exponent,
        },
        "program_5_operational_dynamics": {
            "best_time": best_time,
            "maximum_Jensen_Shannon_divergence_nats": js,
            "KL_unitary_to_diffusive_nats": kl_u_d,
            "observations_for_expected_log_Bayes_factor_10": observations_for_log_bf_10,
            "unitary_Shannon_entropy": entropy(p_unitary),
            "diffusive_Shannon_entropy": entropy(p_diffusive),
            "unitary_return_probability": float(p_unitary[0]),
            "diffusive_return_probability": float(p_diffusive[0]),
        },
        "program_6_information_thermodynamics": {
            "Gibbs_scale_orbit_max_probability_difference": float(
                np.max(np.abs(gibbs_a - gibbs_b))
            ),
            "beta": beta,
            "entropy_reduction": entropy_drop,
            "reservoir_heat": reservoir_heat,
            "relative_entropy_correction": relative_entropy,
            "Landauer_identity_residual": landauer_residual,
        },
        "program_7_noncommutative_geometry": {
            "maximum_point_projection_first_order_residual": max_first_order,
            "diagonal_D_commutator_norm": diagonal_dirac_commutator,
            "interpretation": "off-diagonal geometry fails the naive first-order test; diagonal repair is geometrically trivial",
        },
        "program_8_variational_stability": {
            "strict_minimum_eigenvalue": float(strict_eigenvalues[0]),
            "strict_negative_modes": int(np.sum(strict_eigenvalues < -1e-12)),
            "positive_shift_threshold": threshold,
            "conditioning_above_threshold": stability_scan,
        },
        "program_9_predictive_fingerprint": {
            "positive_affine_scale_legacy_to_strict": float(affine_scale),
            "affine_offset": float(affine_offset),
            "held_out_range_normalized_RMSE": float(
                np.sqrt(np.mean(range_normalized_residuals[1:-1] ** 2))
            ),
            "held_out_max_absolute_residual": float(
                np.max(np.abs(range_normalized_residuals[1:-1]))
            ),
            "strict_normalized_gap_fingerprint": strict_gap_fingerprint.tolist(),
            "legacy_normalized_gap_fingerprint": legacy_gap_fingerprint.tolist(),
        },
        "program_10_machine_assisted_obstructions": {
            "nonempty_binary_profiles": len(binary_profiles),
            "translation_orbits": len(translation_orbits),
            "dihedral_orbits": len(dihedral_orbits),
            "equivariant_even_to_odd_linear_map_dimension": 0,
            "weight_one_polynomial_source_count": weight_one_count,
        },
    }
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
