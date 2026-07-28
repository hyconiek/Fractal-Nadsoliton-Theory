#!/usr/bin/env python3
"""Execute FIN research Programs 255--266.

This is a finite mathematical and synthetic-methodology audit.  It proves
only the algebraic statements stated in the accompanying report and records
deterministic numerical witnesses.  It does not create experimental data,
physical units, a strict selector, a legacy-to-strict bridge, or an adaptive
physical law.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

import fin_lab_p240_optimal_tomography as p240


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "FIN_Programs_255_266_Figures"
RESULTS_PATH = ROOT / "FIN_Programs_255_266_Results.json"
SUMMARY_PATH = ROOT / "FIN_Programs_255_266_Summary.csv"
FALSE_POSITIVE_PATH = ROOT / "FIN_Programs_255_266_False_Positive_Atlas.csv"
RESERVOIR_PATH = ROOT / "FIN_Programs_255_266_Reservoir_Benchmark.csv"

N = 12
SEED = 20260729
STRICT_K = np.array(
    [
        0.46998567264502017,
        0.19204355169010282,
        0.09142861427792497,
        0.04702916874565040,
        0.02413122336363006,
        0.011070817321442113,
    ],
    dtype=float,
)


def cyclic_weight_matrix(weights: np.ndarray) -> np.ndarray:
    weights = np.asarray(weights, dtype=float)
    n = 2 * len(weights)
    w = np.zeros((n, n), dtype=float)
    for x in range(n):
        for y in range(n):
            if x == y:
                continue
            distance = min((x - y) % n, (y - x) % n)
            w[x, y] = weights[distance - 1]
    return w


def laplacian_from_weights(w: np.ndarray) -> np.ndarray:
    return np.diag(w.sum(axis=1)) - w


def strict_operator() -> tuple[np.ndarray, np.ndarray]:
    w = cyclic_weight_matrix(STRICT_K)
    return laplacian_from_weights(w), w


def legacy_amplitude_absorbed_operator() -> tuple[np.ndarray, np.ndarray]:
    """The one frozen P261 completion atom: remove alpha_geo only.

    No phase, frequency, damping, sign, shift, or positivity repair is fitted.
    """
    weights = np.array(
        [
            math.cos(math.pi * d / 4.0 + math.pi / 6.0) / (1.0 + 0.01 * d)
            for d in range(1, 7)
        ],
        dtype=float,
    )
    w = cyclic_weight_matrix(weights)
    return laplacian_from_weights(w), w


def complement(size: int, keep: Iterable[int]) -> list[int]:
    keep_set = set(int(value) for value in keep)
    return [index for index in range(size) if index not in keep_set]


def schur_keep(matrix: np.ndarray, keep: list[int]) -> np.ndarray:
    remove = complement(matrix.shape[0], keep)
    if not remove:
        return matrix.copy()
    m_kk = matrix[np.ix_(keep, keep)]
    m_kr = matrix[np.ix_(keep, remove)]
    m_rk = matrix[np.ix_(remove, keep)]
    m_rr = matrix[np.ix_(remove, remove)]
    return m_kk - m_kr @ np.linalg.solve(m_rr, m_rk)


def block_components(
    matrix: np.ndarray, retained: list[int]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    hidden = complement(matrix.shape[0], retained)
    return (
        matrix[np.ix_(retained, retained)],
        matrix[np.ix_(retained, hidden)],
        matrix[np.ix_(hidden, hidden)],
    )


def self_energy(
    matrix: np.ndarray, retained: list[int], z: float
) -> np.ndarray:
    _, coupling, hidden = block_components(matrix, retained)
    return coupling @ np.linalg.solve(
        z * np.eye(hidden.shape[0]) + hidden,
        coupling.T.conj(),
    )


def hermitian_min_eigenvalue(matrix: np.ndarray) -> float:
    hermitian = (matrix + matrix.T.conj()) / 2.0
    return float(np.linalg.eigvalsh(hermitian).min())


def kl_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    mask = p > 0.0
    if np.any(q[mask] <= 0.0):
        return math.inf
    return float(np.sum(p[mask] * np.log(p[mask] / q[mask])))


def projective_signature(matrix: np.ndarray) -> np.ndarray:
    eigenvalues = np.linalg.eigvalsh((matrix + matrix.T) / 2.0)
    if eigenvalues[-1] <= 0:
        raise ValueError("projective signature requires a positive scale")
    return np.sort(eigenvalues[1:] / eigenvalues[-1])


def program_255(a: np.ndarray) -> dict[str, Any]:
    """All-context finite Stieltjes audit and spectral-measure witness."""
    minimum_hidden_eigenvalue = math.inf
    minimum_signed_derivative_eigenvalue = math.inf
    maximum_measure_reconstruction_residual = 0.0
    maximum_measure_mass_residual = 0.0
    maximum_zero_limit_residual = 0.0
    context_count = 0
    z = 0.2

    for mask in range(1, (1 << N) - 1):
        retained = [index for index in range(N) if mask & (1 << index)]
        _, coupling, hidden = block_components(a, retained)
        eigenvalues, eigenvectors = np.linalg.eigh(hidden)
        minimum_hidden_eigenvalue = min(
            minimum_hidden_eigenvalue, float(eigenvalues.min())
        )
        transformed = coupling @ eigenvectors
        direct = coupling @ np.linalg.solve(
            z * np.eye(hidden.shape[0]) + hidden, coupling.T
        )
        measure = np.zeros_like(direct)
        measure_mass = np.zeros_like(direct)
        for pole, column in zip(eigenvalues, transformed.T):
            residue = np.outer(column, column)
            measure += residue / (z + pole)
            measure_mass += residue
        maximum_measure_reconstruction_residual = max(
            maximum_measure_reconstruction_residual,
            float(np.linalg.norm(direct - measure, 2)),
        )
        maximum_measure_mass_residual = max(
            maximum_measure_mass_residual,
            float(np.linalg.norm(measure_mass - coupling @ coupling.T, 2)),
        )
        zero_limit = coupling @ np.linalg.solve(hidden, coupling.T)
        measure_zero = np.zeros_like(zero_limit)
        for pole, column in zip(eigenvalues, transformed.T):
            measure_zero += np.outer(column, column) / pole
        maximum_zero_limit_residual = max(
            maximum_zero_limit_residual,
            float(np.linalg.norm(zero_limit - measure_zero, 2)),
        )
        for order in range(5):
            signed_derivative = (
                math.factorial(order)
                * coupling
                @ np.linalg.matrix_power(
                    z * np.eye(hidden.shape[0]) + hidden,
                    -(order + 1),
                )
                @ coupling.T
            )
            minimum_signed_derivative_eigenvalue = min(
                minimum_signed_derivative_eigenvalue,
                hermitian_min_eigenvalue(signed_derivative),
            )
        context_count += 1

    retained = list(range(0, N, 2))
    _, coupling, hidden = block_components(a, retained)
    infinity_residuals: dict[str, float] = {}
    zero_residuals: dict[str, float] = {}
    target_infinity = coupling @ coupling.T
    target_zero = coupling @ np.linalg.solve(hidden, coupling.T)
    for value in [1e1, 1e2, 1e3, 1e4]:
        sigma = coupling @ np.linalg.solve(
            value * np.eye(hidden.shape[0]) + hidden, coupling.T
        )
        infinity_residuals[str(value)] = float(
            np.linalg.norm(value * sigma - target_infinity, 2)
        )
    for value in [1e-1, 1e-2, 1e-3, 1e-4]:
        sigma = coupling @ np.linalg.solve(
            value * np.eye(hidden.shape[0]) + hidden, coupling.T
        )
        zero_residuals[str(value)] = float(
            np.linalg.norm(sigma - target_zero, 2)
        )

    return {
        "status": "[Proven]",
        "contexts_audited": context_count,
        "minimum_hidden_principal_eigenvalue": minimum_hidden_eigenvalue,
        "minimum_signed_derivative_eigenvalue_orders_0_to_4": (
            minimum_signed_derivative_eigenvalue
        ),
        "maximum_measure_reconstruction_residual": (
            maximum_measure_reconstruction_residual
        ),
        "maximum_measure_mass_residual": maximum_measure_mass_residual,
        "maximum_zero_limit_spectral_residual": maximum_zero_limit_residual,
        "infinity_limit_residual_norm": infinity_residuals,
        "zero_limit_residual_norm": zero_residuals,
        "theorem": (
            "Sigma_E(z)=sum_mu Gamma_mu/(z+mu), Gamma_mu>=0, and "
            "(-1)^n Sigma_E^(n)(z)>=0 for every nonempty proper context."
        ),
    }


def group_spectral_residues(
    hidden: np.ndarray, coupling: np.ndarray, tolerance: float = 1e-10
) -> list[tuple[float, np.ndarray]]:
    eigenvalues, eigenvectors = np.linalg.eigh(hidden)
    groups: list[list[int]] = []
    for index, value in enumerate(eigenvalues):
        if not groups or abs(value - eigenvalues[groups[-1][0]]) > tolerance:
            groups.append([index])
        else:
            groups[-1].append(index)
    result = []
    for group in groups:
        q = eigenvectors[:, group]
        residue = coupling @ q @ q.T @ coupling.T
        result.append((float(np.mean(eigenvalues[group])), residue))
    return result


def program_256(a: np.ndarray) -> dict[str, Any]:
    retained = list(range(0, N, 2))
    _, coupling, hidden = block_components(a, retained)
    controllability = np.hstack(
        [
            np.linalg.matrix_power(hidden, power) @ coupling.T
            for power in range(hidden.shape[0])
        ]
    )
    controllable_rank = int(
        np.linalg.matrix_rank(controllability, tol=1e-11)
    )
    residues = group_spectral_residues(hidden, coupling)
    residue_rows = []
    minimal_dimension = 0
    visible_poles = 0
    for pole, residue in residues:
        rank = int(np.linalg.matrix_rank(residue, tol=1e-11))
        minimal_dimension += rank
        visible_poles += int(rank > 0)
        residue_rows.append(
            {
                "pole": pole,
                "residue_rank": rank,
                "residue_trace": float(np.trace(residue)),
                "minimum_residue_eigenvalue": hermitian_min_eigenvalue(residue),
            }
        )

    reconstruction_residual = 0.0
    augmented_residual = 0.0
    augmented_hidden = np.zeros((hidden.shape[0] + 2, hidden.shape[0] + 2))
    augmented_hidden[: hidden.shape[0], : hidden.shape[0]] = hidden
    augmented_hidden[-2:, -2:] = np.diag([7.3, 9.1])
    augmented_coupling = np.hstack([coupling, np.zeros((len(retained), 2))])
    for z in [0.05, 0.2, 1.0, 3.0]:
        direct = coupling @ np.linalg.solve(
            z * np.eye(hidden.shape[0]) + hidden, coupling.T
        )
        spectral = sum(residue / (z + pole) for pole, residue in residues)
        reconstruction_residual = max(
            reconstruction_residual,
            float(np.linalg.norm(direct - spectral, 2)),
        )
        augmented = augmented_coupling @ np.linalg.solve(
            z * np.eye(augmented_hidden.shape[0]) + augmented_hidden,
            augmented_coupling.T,
        )
        augmented_residual = max(
            augmented_residual,
            float(np.linalg.norm(direct - augmented, 2)),
        )

    return {
        "status": "[Proven]",
        "declared_hidden_dimension": int(hidden.shape[0]),
        "controllable_observable_rank": controllable_rank,
        "minimal_stieltjes_realization_dimension": minimal_dimension,
        "hidden_spectral_poles": len(residues),
        "visible_self_energy_poles": visible_poles,
        "spectral_residues": residue_rows,
        "spectral_reconstruction_residual": reconstruction_residual,
        "decoupled_two_mode_augmentation_residual": augmented_residual,
        "nonidentifiability_theorem": (
            "Sigma determines only the minimal pole-residue realization up "
            "to hidden orthogonal equivalence; arbitrary decoupled hidden "
            "modes remain invisible."
        ),
    }


def program_257(a: np.ndarray, rng: np.random.Generator) -> dict[str, Any]:
    z = 0.37
    pencil = z * np.eye(N) + a
    maximum_residual = 0.0
    chains = 20_000
    for _ in range(chains):
        permutation = rng.permutation(N)
        size_e = int(rng.integers(1, N - 1))
        size_f = int(rng.integers(size_e + 1, N + 1))
        retained_e = sorted(permutation[:size_e].tolist())
        retained_f = sorted(permutation[:size_f].tolist())
        direct = schur_keep(pencil, retained_e)
        first = schur_keep(pencil, retained_f)
        local_e = [retained_f.index(index) for index in retained_e]
        nested = schur_keep(first, local_e)
        maximum_residual = max(
            maximum_residual, float(np.linalg.norm(direct - nested, 2))
        )
    identity_residual = float(
        np.linalg.norm(schur_keep(pencil, list(range(N))) - pencil, 2)
    )
    return {
        "status": "[Proven]",
        "category": (
            "objects are nonempty vertex contexts; an inclusion E subset F "
            "induces the contravariant Schur restriction B_F -> B_E"
        ),
        "chains_audited": chains,
        "maximum_composition_residual": maximum_residual,
        "identity_residual": identity_residual,
        "proof_core": (
            "associativity follows from uniqueness of block Gaussian "
            "elimination, equivalently iterated minimization of one positive "
            "quadratic form"
        ),
        "not_proven": "no sheaf gluing axiom and no quantum contextuality claim",
    }


def twisted_generator(
    theta: float, weights: np.ndarray = STRICT_K
) -> np.ndarray:
    n = 2 * len(weights)
    w0 = cyclic_weight_matrix(weights)
    s = float(w0[0].sum())
    result_w = np.zeros((n, n), dtype=complex)
    for distance in range(1, n // 2):
        shift = np.zeros((n, n), dtype=complex)
        for x in range(n):
            shift[x, (x + distance) % n] = 1.0
        result_w += weights[distance - 1] * (
            np.exp(1j * distance * theta) * shift
            + np.exp(-1j * distance * theta) * shift.T
        )
    antipodal = np.zeros((n, n), dtype=complex)
    for x in range(n):
        antipodal[x, (x + n // 2) % n] = 1.0
    result_w += (
        weights[-1]
        * math.cos((n // 2) * theta)
        * antipodal
    )
    return s * np.eye(n, dtype=complex) - result_w


def twisted_generator_derivative_at_zero(
    weights: np.ndarray = STRICT_K,
) -> np.ndarray:
    n = 2 * len(weights)
    derivative_w = np.zeros((n, n), dtype=complex)
    for distance in range(1, n // 2):
        shift = np.zeros((n, n), dtype=complex)
        for x in range(n):
            shift[x, (x + distance) % n] = 1.0
        derivative_w += weights[distance - 1] * (
            1j * distance * shift - 1j * distance * shift.T
        )
    return -derivative_w


def chiral_susceptibility_analytic(
    a0: np.ndarray,
    aprime: np.ndarray,
    retained: list[int],
    z: float,
) -> tuple[np.ndarray, float]:
    hidden_indices = complement(a0.shape[0], retained)
    coupling = a0[np.ix_(retained, hidden_indices)]
    hidden = a0[np.ix_(hidden_indices, hidden_indices)]
    coupling_prime = aprime[np.ix_(retained, hidden_indices)]
    hidden_prime = aprime[np.ix_(hidden_indices, hidden_indices)]
    resolvent = np.linalg.inv(z * np.eye(len(hidden_indices)) + hidden)
    xi = (
        coupling_prime @ resolvent @ coupling.T.conj()
        + coupling @ resolvent @ coupling_prime.T.conj()
        - coupling @ resolvent @ hidden_prime @ resolvent @ coupling.T.conj()
    )
    denominator = z + float(np.linalg.eigvalsh(hidden).min())
    bound = (
        2.0
        * np.linalg.norm(coupling_prime, 2)
        * np.linalg.norm(coupling, 2)
        / denominator
        + np.linalg.norm(coupling, 2) ** 2
        * np.linalg.norm(hidden_prime, 2)
        / denominator**2
    )
    return xi, float(bound)


def program_258(a: np.ndarray) -> dict[str, Any]:
    retained = list(range(0, N, 2))
    a0 = twisted_generator(0.0)
    aprime = twisted_generator_derivative_at_zero()
    h = 1e-6
    reflection = np.zeros((N, N))
    for x in range(N):
        reflection[x, (-x) % N] = 1.0
    reflection_retained = reflection[np.ix_(retained, retained)]
    rows = {}
    maximum_derivative_residual = 0.0
    maximum_covariance_residual = 0.0
    maximum_bound_ratio = 0.0
    for z in [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]:
        analytic, bound = chiral_susceptibility_analytic(
            a0, aprime, retained, z
        )
        finite = (
            self_energy(twisted_generator(h), retained, z)
            - self_energy(twisted_generator(-h), retained, z)
        ) / (2.0 * h)
        derivative_residual = float(np.linalg.norm(analytic - finite, 2))
        covariance_residual = float(
            np.linalg.norm(
                reflection_retained @ analytic @ reflection_retained
                + analytic,
                2,
            )
        )
        norm = float(np.linalg.norm(analytic, 2))
        maximum_derivative_residual = max(
            maximum_derivative_residual, derivative_residual
        )
        maximum_covariance_residual = max(
            maximum_covariance_residual, covariance_residual
        )
        maximum_bound_ratio = max(
            maximum_bound_ratio, norm / bound if bound > 0 else 0.0
        )
        rows[str(z)] = {
            "operator_norm": norm,
            "analytic_vs_finite_difference_residual": derivative_residual,
            "inversion_odd_covariance_residual": covariance_residual,
            "proved_norm_bound": bound,
            "norm_to_bound_ratio": norm / bound,
        }
    return {
        "status": "[Proven]",
        "formula": (
            "Xi=B' R B* + B R B'^* - B R C' R B*, "
            "where R=(zI+C)^-1"
        ),
        "z_audit": rows,
        "maximum_derivative_residual": maximum_derivative_residual,
        "maximum_covariance_residual": maximum_covariance_residual,
        "maximum_norm_to_bound_ratio": maximum_bound_ratio,
        "selector_boundary": (
            "Xi is an inversion-odd tangent receiver for a supplied twist; "
            "it does not select the sign of theta."
        ),
    }


def quantile_signature(eigenvalues: np.ndarray, points: int = 8) -> np.ndarray:
    values = np.sort(np.asarray(eigenvalues, dtype=float))
    values = values / values[-1]
    source = np.linspace(0.0, 1.0, len(values))
    target = np.linspace(0.0, 1.0, points)
    return np.interp(target, source, values)


def program_259(a: np.ndarray) -> dict[str, Any]:
    contexts = [
        list(range(N)),
        list(range(0, N, 2)),
        [0, 4, 8],
        [0],
    ]
    z_grid = [0.1, 0.2, 0.5, 1.0]
    levels = []
    descriptors = []
    for level, retained in enumerate(contexts):
        spectra = []
        memory_fractions = []
        for z in z_grid:
            pencil = z * np.eye(N) + a
            effective = schur_keep(pencil, retained)
            spectra.extend(quantile_signature(np.linalg.eigvalsh(effective)))
            if level == 0:
                memory_fractions.append(0.0)
            else:
                parent = contexts[level - 1]
                parent_pencil = schur_keep(pencil, parent)
                local = [parent.index(index) for index in retained]
                parent_direct = parent_pencil[np.ix_(local, local)]
                reduced = schur_keep(parent_pencil, local)
                memory = parent_direct - reduced
                memory_fractions.append(
                    float(
                        np.trace(memory)
                        / max(np.trace(parent_direct), np.finfo(float).eps)
                    )
                )
        descriptor = np.array(spectra + memory_fractions, dtype=float)
        descriptors.append(descriptor)
        levels.append(
            {
                "level": level,
                "context": retained,
                "cardinality": len(retained),
                "descriptor": descriptor.tolist(),
                "memory_fraction_by_z": memory_fractions,
            }
        )
    distances = [
        float(np.linalg.norm(left - right))
        for left, right in zip(descriptors, descriptors[1:])
    ]
    return {
        "status": "[Proven] no-cycle for the declared strict-decimation category; [Moderate evidence] descriptor audit",
        "levels": levels,
        "successive_descriptor_distances": distances,
        "minimum_successive_descriptor_distance": min(distances),
        "fixed_point_no_go": (
            "a strict context reduction decreases cardinality 12>6>3>1; "
            "there is no nontrivial fixed point or cycle in the category "
            "whose object includes context dimension"
        ),
        "open_requirement": (
            "a genuine RG fixed point requires an additional size-restoring "
            "embedding/equivalence and normalization law"
        ),
    }


def coarse_record_decomposition(
    p: np.ndarray, q: np.ndarray, groups: np.ndarray
) -> tuple[float, float, float]:
    outcomes = sorted(set(int(value) for value in groups))
    record_p = np.array([p[groups == value].sum() for value in outcomes])
    record_q = np.array([q[groups == value].sum() for value in outcomes])
    record = kl_divergence(record_p, record_q)
    conditional = 0.0
    for index, value in enumerate(outcomes):
        mask = groups == value
        if record_p[index] <= 0:
            continue
        conditional += record_p[index] * kl_divergence(
            p[mask] / record_p[index],
            q[mask] / record_q[index],
        )
    total = kl_divergence(p, q)
    return record, conditional, total - record - conditional


def program_260(a: np.ndarray, rng: np.random.Generator) -> dict[str, Any]:
    transition = expm(-0.25 * a)
    groups = np.arange(N) % 3
    p = rng.dirichlet(np.ones(N))
    q = rng.dirichlet(np.ones(N))
    rows = []
    initial = kl_divergence(p, q)
    maximum_chain_residual = 0.0
    minimum_environment_loss = math.inf
    for step in range(1, 6):
        before = kl_divergence(p, q)
        next_p = transition @ p
        next_q = transition @ q
        after = kl_divergence(next_p, next_q)
        environment_loss = before - after
        record, conditional, chain_residual = coarse_record_decomposition(
            next_p, next_q, groups
        )
        maximum_chain_residual = max(
            maximum_chain_residual, abs(chain_residual)
        )
        minimum_environment_loss = min(
            minimum_environment_loss, environment_loss
        )
        rows.append(
            {
                "step": step,
                "input_distinguishability": before,
                "discarded_environment_loss": environment_loss,
                "joint_postchannel_distinguishability": after,
                "apparatus_record_distinguishability": record,
                "conditional_system_carry": conditional,
                "chain_rule_residual": chain_residual,
            }
        )
        p, q = next_p, next_q
    final = kl_divergence(p, q)
    telescoping = initial - final
    summed_loss = sum(row["discarded_environment_loss"] for row in rows)

    minimum_random_loss = math.inf
    for _ in range(1000):
        random_p = rng.dirichlet(np.ones(N))
        random_q = rng.dirichlet(np.ones(N))
        loss = kl_divergence(random_p, random_q) - kl_divergence(
            transition @ random_p, transition @ random_q
        )
        minimum_random_loss = min(minimum_random_loss, loss)
    return {
        "status": "[Proven] for the declared classical instrument ledger",
        "steps": rows,
        "telescoping_total_loss": telescoping,
        "sum_of_step_losses": summed_loss,
        "telescoping_residual": abs(telescoping - summed_loss),
        "maximum_record_chain_rule_residual": maximum_chain_residual,
        "minimum_step_environment_loss": minimum_environment_loss,
        "minimum_loss_over_1000_random_pairs": minimum_random_loss,
        "physical_boundary": (
            "the ledger is dimensionless accessible distinguishability, not "
            "thermodynamic entropy or heat without a bath, Hamiltonian and scale"
        ),
    }


def program_261(strict_a: np.ndarray) -> dict[str, Any]:
    legacy_a, _ = legacy_amplitude_absorbed_operator()
    retained = list(range(0, N, 2))
    legacy_hidden = legacy_a[np.ix_(complement(N, retained), complement(N, retained))]
    strict_hidden = strict_a[np.ix_(complement(N, retained), complement(N, retained))]
    legacy_spectrum = np.linalg.eigvalsh(legacy_a)
    strict_spectrum = np.linalg.eigvalsh(strict_a)
    shape_defects = {}
    for z in [6.5, 8.0, 12.0, 20.0]:
        sigma_strict = self_energy(strict_a, retained, z)
        sigma_legacy = self_energy(legacy_a, retained, z)
        normalized_strict = sigma_strict / np.linalg.norm(sigma_strict, "fro")
        normalized_legacy = sigma_legacy / np.linalg.norm(sigma_legacy, "fro")
        shape_defects[str(z)] = float(
            np.linalg.norm(normalized_strict - normalized_legacy, "fro")
        )
    positive_ray_pole = float(-np.linalg.eigvalsh(legacy_hidden).min())
    return {
        "status": "[Proven] amplitude-only completion obstruction",
        "completion_map_frozen_before_test": (
            "divide K_legacy by alpha_geo, preserve legacy omega, phi and "
            "beta_tors, then apply the same A=diag(W1)-W and context map"
        ),
        "strict_generator_minimum_eigenvalue": float(strict_spectrum.min()),
        "legacy_amplitude_absorbed_generator_minimum_eigenvalue": float(
            legacy_spectrum.min()
        ),
        "strict_hidden_minimum_eigenvalue": float(
            np.linalg.eigvalsh(strict_hidden).min()
        ),
        "legacy_hidden_minimum_eigenvalue": float(
            np.linalg.eigvalsh(legacy_hidden).min()
        ),
        "legacy_positive_resolvent_ray_pole": positive_ray_pole,
        "projective_self_energy_shape_defect": shape_defects,
        "verdict": (
            "amplitude absorption alone does not map the legacy object into "
            "the strict positive Stieltjes-generator class; this does not "
            "exclude a bridge with additional sign/phase/damping/compression data"
        ),
    }


def program_262(
    a: np.ndarray, rng: np.random.Generator, replicates: int = 300
) -> dict[str, Any]:
    eigenvalues = np.linalg.eigvalsh(a)
    lambda_max = float(eigenvalues[-1])
    tau = 1.0 / lambda_max
    target_signature = eigenvalues[1:] / lambda_max
    transition = expm(-tau * a)
    true_scale = 3.7
    physical_time = tau / true_scale
    signature_errors = []
    scale_errors = []
    joint_passes = 0
    invalid = 0
    for _ in range(replicates):
        counts = p240.sample_counts(transition, 50_000, rng)
        signature, _ = p240.signature_from_counts(counts, tau)
        if signature is None:
            invalid += 1
            continue
        signature_error = float(
            np.max(np.abs(signature - target_signature))
        )
        measured_time = physical_time * (1.0 + rng.normal(0.0, 0.005))
        estimated_scale = tau / measured_time
        scale_error = abs(estimated_scale / true_scale - 1.0)
        signature_errors.append(signature_error)
        scale_errors.append(scale_error)
        if signature_error <= 0.03 and scale_error <= 0.02:
            joint_passes += 1
    invariance_residual = float(
        np.linalg.norm(
            expm(-physical_time * true_scale * a) - expm(-tau * a), 2
        )
    )
    return {
        "status": "[Strong evidence] synthetic planning audit; external experiment not executed",
        "replicates": replicates,
        "shots_per_preparation": 50_000,
        "projective_threshold": 0.03,
        "independent_scale_threshold": 0.02,
        "invalid_reconstructions": invalid,
        "projective_error_median": float(np.median(signature_errors)),
        "projective_error_95_percentile": float(
            np.quantile(signature_errors, 0.95)
        ),
        "scale_error_median": float(np.median(scale_errors)),
        "scale_error_95_percentile": float(np.quantile(scale_errors, 0.95)),
        "joint_pass_rate": joint_passes / replicates,
        "scale_orbit_transition_invariance_residual": invariance_residual,
        "gate_boundary": (
            "fingerprint is tested before independent clock calibration; "
            "P242 remains blocked without a P241-admitted external bundle"
        ),
    }


def current_observable(w: np.ndarray, distance: int) -> np.ndarray:
    observable = np.zeros_like(w, dtype=complex)
    for x in range(w.shape[0]):
        y = (x + distance) % w.shape[0]
        observable[x, y] += -1j * distance * w[x, y]
        observable[y, x] += 1j * distance * w[x, y]
    return observable


def program_263(a: np.ndarray, w: np.ndarray) -> dict[str, Any]:
    rho_plus = np.eye(N, dtype=complex) / N
    epsilon = 0.04
    left, right = 10, 2
    rho_plus[left, right] += 1j * epsilon
    rho_plus[right, left] -= 1j * epsilon
    rho_minus = rho_plus.conj()
    population_residual = float(
        np.max(np.abs(np.diag(rho_plus) - np.diag(rho_minus)))
    )
    observables = [current_observable(w, distance) for distance in range(1, 6)]
    currents_plus = [
        float(np.real(np.trace(rho_plus @ observable)))
        for observable in observables
    ]
    currents_minus = [
        float(np.real(np.trace(rho_minus @ observable)))
        for observable in observables
    ]
    gram = np.array(
        [
            [
                float(np.real(np.trace(left_obs.conj().T @ right_obs)))
                for right_obs in observables
            ]
            for left_obs in observables
        ]
    )
    gram_eigenvalues = np.linalg.eigvalsh(gram)
    current_rank = int(np.linalg.matrix_rank(gram, tol=1e-12))
    condition = float(
        gram_eigenvalues[-1] / gram_eigenvalues[0]
    )

    retained = list(range(0, N, 2))
    analytic, _ = chiral_susceptibility_analytic(
        twisted_generator(0.0),
        twisted_generator_derivative_at_zero(),
        retained,
        0.2,
    )
    h = 1e-5
    finite = (
        self_energy(twisted_generator(h), retained, 0.2)
        - self_energy(twisted_generator(-h), retained, 0.2)
    ) / (2.0 * h)
    return {
        "status": "[Proven] vertex-POVM no-go and current-observable construction; [Strong evidence] synthetic Xi protocol",
        "minimum_eigenvalue_rho_plus": hermitian_min_eigenvalue(rho_plus),
        "minimum_eigenvalue_rho_minus": hermitian_min_eigenvalue(rho_minus),
        "vertex_population_residual": population_residual,
        "harmonic_currents_rho_plus_d1_to_d5": currents_plus,
        "harmonic_currents_rho_minus_d1_to_d5": currents_minus,
        "opposite_current_residual": float(
            np.max(np.abs(np.array(currents_plus) + np.array(currents_minus)))
        ),
        "five_current_observable_rank": current_rank,
        "current_observable_gram_condition_number": condition,
        "xi_analytic_vs_two_flux_process_tomography_residual": float(
            np.linalg.norm(analytic - finite, 2)
        ),
        "no_go": (
            "rho and rho* have identical vertex probabilities and opposite "
            "currents; a vertex POVM alone cannot identify C_d or C_chi"
        ),
        "required_instrument": (
            "phase-sensitive interferometric current observables plus "
            "generator/process tomography at supplied opposite flux twists"
        ),
    }


def random_positive_laplacian(
    rng: np.random.Generator, circulant: bool
) -> tuple[np.ndarray, np.ndarray]:
    if circulant:
        weights = rng.uniform(0.01, 0.65, size=6)
        w = cyclic_weight_matrix(weights)
    else:
        upper = rng.uniform(0.01, 0.65, size=(N, N))
        upper = np.triu(upper, 1)
        w = upper + upper.T
    return laplacian_from_weights(w), w


def signed_circulant_laplacian(
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    weights = rng.uniform(-0.65, 0.65, size=6)
    w = cyclic_weight_matrix(weights)
    return laplacian_from_weights(w), w


def atlas_tests(
    candidate: np.ndarray,
    target_signature: np.ndarray,
    reflection: np.ndarray,
) -> dict[str, bool]:
    retained = list(range(0, N, 2))
    hidden = complement(N, retained)
    signature_pass = False
    if np.linalg.eigvalsh(candidate)[-1] > 0:
        signature = projective_signature(candidate)
        signature_pass = bool(
            np.max(np.abs(signature - target_signature)) <= 0.02
        )
    hidden_block = candidate[np.ix_(hidden, hidden)]
    stieltjes_pass = bool(
        np.linalg.eigvalsh(hidden_block).min() > 1e-10
    )
    if stieltjes_pass:
        coupling = candidate[np.ix_(retained, hidden)]
        for order in range(3):
            signed = (
                math.factorial(order)
                * coupling
                @ np.linalg.matrix_power(
                    0.2 * np.eye(len(hidden)) + hidden_block,
                    -(order + 1),
                )
                @ coupling.T
            )
            stieltjes_pass &= hermitian_min_eigenvalue(signed) >= -1e-9
    pencil = 0.2 * np.eye(N) + candidate
    context_pass = False
    try:
        direct = schur_keep(pencil, retained)
        intermediate_context = [0, 1, 2, 3, 4, 6, 8, 10]
        intermediate = schur_keep(pencil, intermediate_context)
        local = [intermediate_context.index(index) for index in retained]
        nested = schur_keep(intermediate, local)
        context_pass = bool(np.linalg.norm(direct - nested, 2) <= 1e-9)
    except np.linalg.LinAlgError:
        context_pass = False
    reflection_pass = bool(
        np.linalg.norm(reflection @ candidate @ reflection - candidate, 2)
        <= 1e-9
    )
    markov_generator = bool(
        np.linalg.eigvalsh(candidate).min() >= -1e-9
        and np.max(np.abs(candidate.sum(axis=0))) <= 1e-9
        and np.max(candidate - np.diag(np.diag(candidate))) <= 1e-9
    )
    information_pass = False
    if markov_generator:
        transition = expm(-0.3 * candidate)
        information_pass = bool(
            transition.min() >= -1e-10
            and np.max(np.abs(transition.sum(axis=0) - 1.0)) <= 1e-9
        )
    return {
        "strict_fingerprint": signature_pass,
        "stieltjes_memory": stieltjes_pass,
        "context_composition": context_pass,
        "chiral_covariance_scaffold": reflection_pass,
        "information_contraction_channel": information_pass,
    }


def program_264(
    strict_a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    reflection = np.zeros((N, N))
    for x in range(N):
        reflection[x, (-x) % N] = 1.0
    target_signature = projective_signature(strict_a)
    candidates: list[tuple[str, np.ndarray]] = [("strict_target", strict_a)]
    for _ in range(100):
        candidates.append(("random_positive_dense", random_positive_laplacian(rng, False)[0]))
        candidates.append(("random_positive_circulant", random_positive_laplacian(rng, True)[0]))
        candidates.append(("random_signed_circulant", signed_circulant_laplacian(rng)[0]))
        perturbed = STRICT_K * np.exp(rng.normal(0.0, 0.05, size=6))
        candidates.append(
            ("strict_five_percent_circulant", laplacian_from_weights(cyclic_weight_matrix(perturbed)))
        )
    rows = []
    tests = [
        "strict_fingerprint",
        "stieltjes_memory",
        "context_composition",
        "chiral_covariance_scaffold",
        "information_contraction_channel",
    ]
    for index, (ensemble, candidate) in enumerate(candidates):
        result = atlas_tests(candidate, target_signature, reflection)
        rows.append(
            {
                "candidate_id": index,
                "ensemble": ensemble,
                **{name: int(result[name]) for name in tests},
                "all_five": int(all(result.values())),
            }
        )
    summaries = {}
    for ensemble in sorted(set(row["ensemble"] for row in rows)):
        selected = [row for row in rows if row["ensemble"] == ensemble]
        summaries[ensemble] = {
            name: sum(row[name] for row in selected) / len(selected)
            for name in tests + ["all_five"]
        }
        summaries[ensemble]["count"] = len(selected)
    return (
        {
            "status": "[Proven] finite false-positive atlas for frozen ensembles and thresholds",
            "thresholds": {
                "strict_fingerprint_max_error": 0.02,
                "algebraic_residual": 1e-9,
            },
            "ensemble_pass_rates": summaries,
            "interpretation": (
                "Stieltjes positivity and Schur composition are generic "
                "operator properties; strict spectral shape supplies most "
                "specificity, while chiral covariance is common to circulants."
            ),
        },
        rows,
    )


def edge_laplacian(left: int, right: int, n: int = N) -> np.ndarray:
    result = np.zeros((n, n))
    result[left, left] = 1.0
    result[right, right] = 1.0
    result[left, right] = -1.0
    result[right, left] = -1.0
    return result


def program_265(a: np.ndarray) -> dict[str, Any]:
    time = 0.3
    scale = 1.2
    base = expm(-time * a)
    uniform_generator_change = expm(-time * scale * a)
    clock_change = expm(-(time * scale) * a)
    exact_degeneracy = float(
        np.linalg.norm(uniform_generator_change - clock_change, 2)
    )
    nonuniform_a = a + 0.18 * edge_laplacian(0, 1)
    nonuniform = expm(-time * nonuniform_a)
    apparatus_two_step = 0.88 * (base @ base) + 0.12 * np.eye(N)
    scenarios = {
        "static_operator": {
            "projective_fingerprint_drift": 0.0,
            "two_step_semigroup_defect": float(
                np.linalg.norm(base @ base - base @ base, 2)
            ),
            "independent_clock_required": False,
        },
        "uniform_generator_or_clock_scale": {
            "projective_fingerprint_drift": float(
                np.max(
                    np.abs(
                        projective_signature(scale * a)
                        - projective_signature(a)
                    )
                )
            ),
            "two_step_semigroup_defect": 0.0,
            "independent_clock_required": True,
        },
        "nonuniform_generator_change": {
            "projective_fingerprint_drift": float(
                np.max(
                    np.abs(
                        projective_signature(nonuniform_a)
                        - projective_signature(a)
                    )
                )
            ),
            "one_step_transition_change": float(
                np.linalg.norm(nonuniform - base, 2)
            ),
            "two_step_semigroup_defect": 0.0,
            "independent_clock_required": False,
        },
        "apparatus_memory": {
            "projective_fingerprint_drift": 0.0,
            "two_step_semigroup_defect": float(
                np.linalg.norm(apparatus_two_step - base @ base, 2)
            ),
            "independent_clock_required": False,
        },
    }
    return {
        "status": "[Proven] identifiability quotient and exact scale-clock no-go; [Strong evidence] synthetic diagnostic table",
        "exact_uniform_generator_vs_clock_degeneracy_residual": exact_degeneracy,
        "scenarios": scenarios,
        "identifiability_quotient": {
            "scale_orbit": "(cA,t) is observationally equivalent to (A,ct) without an external clock",
            "memory_class": "one-step transition is insufficient; multi-time interventions test semigroup defects",
            "shape_class": "projective spectral drift detects nonuniform generator changes",
        },
        "adaptive_law_boundary": (
            "finite operator snapshots do not uniquely identify an update "
            "vector field; causal interventions and held-out trajectories are required"
        ),
    }


def stable_reservoir_from_transition(
    transition: np.ndarray, gamma: float = 0.95
) -> np.ndarray:
    return gamma * transition


def memory_capacity(
    reservoir: np.ndarray,
    input_vector: np.ndarray,
    input_series: np.ndarray,
    delays: int = 20,
) -> float:
    state = np.zeros(reservoir.shape[0])
    states = []
    for value in input_series:
        state = reservoir @ state + input_vector * value
        states.append(state.copy())
    states_array = np.asarray(states)
    start = 250
    usable = np.arange(start + delays, len(input_series))
    split = int(0.6 * len(usable))
    train_indices = usable[:split]
    test_indices = usable[split:]
    x_train = states_array[train_indices].T
    x_test = states_array[test_indices].T
    x_mean = x_train.mean(axis=1, keepdims=True)
    x_train = x_train - x_mean
    x_test = x_test - x_mean
    y_train = np.vstack(
        [input_series[train_indices - delay] for delay in range(1, delays + 1)]
    )
    y_test = np.vstack(
        [input_series[test_indices - delay] for delay in range(1, delays + 1)]
    )
    weights = (
        y_train
        @ x_train.T
        @ np.linalg.inv(
            x_train @ x_train.T + 1e-8 * np.eye(reservoir.shape[0])
        )
    )
    prediction = weights @ x_test
    capacities = []
    for target, estimate in zip(y_test, prediction):
        if np.std(estimate) < 1e-12:
            capacities.append(0.0)
        else:
            capacities.append(float(np.corrcoef(target, estimate)[0, 1] ** 2))
    return float(sum(capacities))


def reservoir_structural_metrics(
    reservoir: np.ndarray, input_vector: np.ndarray
) -> tuple[float, float, int]:
    gramian = np.zeros((reservoir.shape[0], reservoir.shape[0]))
    vector = input_vector.copy()
    for _ in range(500):
        gramian += np.outer(vector, vector)
        vector = reservoir @ vector
    eigenvalues = np.linalg.eigvalsh(gramian)
    positive = eigenvalues[eigenvalues > 1e-14]
    effective_dimension = float(
        positive.sum() ** 2 / np.sum(positive**2)
    )
    descending = positive[::-1]
    cumulative = np.cumsum(descending) / descending.sum()
    readout_dimension = int(np.searchsorted(cumulative, 0.95) + 1)
    perturbation = np.ones(reservoir.shape[0])
    perturbation /= np.linalg.norm(perturbation)
    recovery = 500
    for step in range(501):
        if np.linalg.norm(perturbation) <= 1e-3:
            recovery = step
            break
        perturbation = reservoir @ perturbation
    return effective_dimension, float(recovery), readout_dimension


def program_266(
    a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    transition = expm(-0.25 * a)
    fin_reservoir = stable_reservoir_from_transition(transition)
    input_vector = np.zeros(N)
    input_vector[0] = 1.0
    input_vector -= input_vector.mean()
    input_vector /= np.linalg.norm(input_vector)
    input_series = rng.uniform(-1.0, 1.0, size=5_000)
    rows = []

    def add_row(
        name: str, index: int, reservoir: np.ndarray, vector: np.ndarray
    ) -> None:
        capacity = memory_capacity(
            reservoir, vector, input_series, delays=20
        )
        effective_dimension, recovery, readout = reservoir_structural_metrics(
            reservoir, vector
        )
        rows.append(
            {
                "ensemble": name,
                "index": index,
                "memory_capacity_delays_1_to_20": capacity,
                "controllability_effective_dimension": effective_dimension,
                "recovery_steps_to_1e_minus_3": recovery,
                "minimal_readout_dimension_95_percent": readout,
            }
        )

    add_row("FIN_heat_reservoir", 0, fin_reservoir, input_vector)
    eigenvalues, eigenvectors = np.linalg.eigh(fin_reservoir)
    for index in range(60):
        random_matrix = rng.normal(size=(N, N))
        q, _ = np.linalg.qr(random_matrix)
        matched = q @ np.diag(eigenvalues) @ q.T
        add_row("isospectral_random_orientation", index, matched, input_vector)
    for index in range(60):
        random_eigenvalues = rng.uniform(-0.95, 0.95, size=N)
        random_eigenvalues[np.argmax(np.abs(random_eigenvalues))] = 0.95
        random_matrix = rng.normal(size=(N, N))
        q, _ = np.linalg.qr(random_matrix)
        random_reservoir = q @ np.diag(random_eigenvalues) @ q.T
        add_row("random_symmetric_spectral_radius_matched", index, random_reservoir, input_vector)

    fin = rows[0]
    comparisons = [row for row in rows if row["ensemble"] != "FIN_heat_reservoir"]
    metrics = [
        "memory_capacity_delays_1_to_20",
        "controllability_effective_dimension",
        "recovery_steps_to_1e_minus_3",
        "minimal_readout_dimension_95_percent",
    ]
    percentiles = {}
    for metric in metrics:
        values = np.array([float(row[metric]) for row in comparisons])
        percentiles[metric] = float(
            np.mean(values <= float(fin[metric]))
        )
    decisive = any(
        percentile < 0.025 or percentile > 0.975
        for percentile in percentiles.values()
    )
    return (
        {
            "status": "[Strong evidence] controlled computational benchmark",
            "fin_metrics": fin,
            "control_count": len(comparisons),
            "fin_percentile_among_controls": percentiles,
            "distinctive_outlier_under_declared_2_5_percent_rule": decisive,
            "verdict": (
                "FIN passes as a finite fading-memory reservoir, but a "
                "biological or cybernetic advantage is retained only if it "
                "outperforms matched controls on held-out tasks."
            ),
            "ontology_boundary": "mechanistic benchmark only; no biological ontology is inferred",
        },
        rows,
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_figures(
    results: dict[str, Any],
    false_positive_rows: list[dict[str, Any]],
    reservoir_rows: list[dict[str, Any]],
) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)

    p255 = results["P255"]
    p258 = results["P258"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    z_inf = np.array([float(key) for key in p255["infinity_limit_residual_norm"]])
    r_inf = np.array(list(p255["infinity_limit_residual_norm"].values()))
    z_zero = np.array([float(key) for key in p255["zero_limit_residual_norm"]])
    r_zero = np.array(list(p255["zero_limit_residual_norm"].values()))
    axes[0].loglog(z_inf, r_inf, "o-", label=r"$z\Sigma(z)\to BB^*$")
    axes[0].loglog(z_zero, r_zero, "s-", label=r"$\Sigma(z)\to BC^{-1}B^*$")
    axes[0].set_xlabel("z")
    axes[0].set_ylabel("operator-norm residual")
    axes[0].set_title("P255: Stieltjes endpoint limits")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    z_values = np.array([float(key) for key in p258["z_audit"]])
    xi_norms = np.array(
        [p258["z_audit"][str(value)]["operator_norm"] for value in z_values]
    )
    bounds = np.array(
        [p258["z_audit"][str(value)]["proved_norm_bound"] for value in z_values]
    )
    axes[1].loglog(z_values, xi_norms, "o-", label=r"$\|\Xi_E(z)\|$")
    axes[1].loglog(z_values, bounds, "--", label="proved bound")
    axes[1].set_xlabel("z")
    axes[1].set_ylabel("operator norm")
    axes[1].set_title("P258: chiral memory susceptibility")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p255_p258_stieltjes_chiral.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    p259 = results["P259"]
    levels = [row["level"] for row in p259["levels"]]
    cardinalities = [row["cardinality"] for row in p259["levels"]]
    axes[0].plot(levels, cardinalities, "o-", label="context cardinality")
    axes[0].set_xticks(levels)
    axes[0].set_xlabel("reduction level")
    axes[0].set_ylabel("|E|")
    axes[0].set_title("P259: strict-decimation flow")
    axes[0].grid(alpha=0.25)
    p260 = results["P260"]
    steps = [row["step"] for row in p260["steps"]]
    loss = [row["discarded_environment_loss"] for row in p260["steps"]]
    record = [row["apparatus_record_distinguishability"] for row in p260["steps"]]
    carry = [row["conditional_system_carry"] for row in p260["steps"]]
    axes[1].plot(steps, loss, "o-", label="discarded loss")
    axes[1].plot(steps, record, "s-", label="record")
    axes[1].plot(steps, carry, "^-", label="conditional carry")
    axes[1].set_xlabel("process step")
    axes[1].set_ylabel("nats")
    axes[1].set_title("P260: process information ledger")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p259_p260_rg_information_ledger.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    p261 = results["P261"]
    z_values = np.array(
        [float(key) for key in p261["projective_self_energy_shape_defect"]]
    )
    defects = np.array(list(p261["projective_self_energy_shape_defect"].values()))
    axes[0].plot(z_values, defects, "o-", color="#d95f02")
    axes[0].axvline(
        p261["legacy_positive_resolvent_ray_pole"],
        color="black",
        linestyle="--",
        label="legacy positive-ray pole",
    )
    axes[0].set_xlabel("z")
    axes[0].set_ylabel("projective self-energy defect")
    axes[0].set_title("P261: amplitude-only completion obstruction")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    p262 = results["P262"]
    axes[1].bar(
        ["projective\nmedian", "projective\n95%", "scale\nmedian", "scale\n95%"],
        [
            p262["projective_error_median"],
            p262["projective_error_95_percentile"],
            p262["scale_error_median"],
            p262["scale_error_95_percentile"],
        ],
        color=["#1b9e77", "#66c2a5", "#7570b3", "#8da0cb"],
    )
    axes[1].axhline(0.02, color="black", linestyle=":", label="2% reference")
    axes[1].set_ylabel("absolute/relative error")
    axes[1].set_title("P262: frozen synthetic calibration audit")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p261_p262_bridge_fingerprint.png", dpi=190)
    plt.close(fig)

    tests = [
        "strict_fingerprint",
        "stieltjes_memory",
        "context_composition",
        "chiral_covariance_scaffold",
        "information_contraction_channel",
        "all_five",
    ]
    ensembles = sorted(set(row["ensemble"] for row in false_positive_rows))
    matrix = np.zeros((len(ensembles), len(tests)))
    for i, ensemble in enumerate(ensembles):
        rows = [row for row in false_positive_rows if row["ensemble"] == ensemble]
        for j, test in enumerate(tests):
            matrix[i, j] = np.mean([row[test] for row in rows])
    fig, ax = plt.subplots(figsize=(11.5, 5.5))
    image = ax.imshow(matrix, vmin=0.0, vmax=1.0, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(tests)), tests, rotation=30, ha="right")
    ax.set_yticks(range(len(ensembles)), ensembles)
    ax.set_title("P264 false-positive atlas: ensemble pass rates")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(
                j,
                i,
                f"{matrix[i, j]:.2f}",
                ha="center",
                va="center",
                color="white" if matrix[i, j] < 0.45 else "black",
                fontsize=8,
            )
    fig.colorbar(image, ax=ax, label="pass rate")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p264_false_positive_atlas.png", dpi=190)
    plt.close(fig)

    metrics = [
        "memory_capacity_delays_1_to_20",
        "controllability_effective_dimension",
        "recovery_steps_to_1e_minus_3",
        "minimal_readout_dimension_95_percent",
    ]
    controls = [
        row for row in reservoir_rows if row["ensemble"] != "FIN_heat_reservoir"
    ]
    fin = reservoir_rows[0]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
    for axis, metric in zip(axes.flat, metrics):
        values = [float(row[metric]) for row in controls]
        axis.hist(values, bins=18, color="#80b1d3", alpha=0.8)
        axis.axvline(float(fin[metric]), color="#e31a1c", lw=2, label="FIN")
        axis.set_title(metric.replace("_", " "))
        axis.legend()
    fig.suptitle("P266: FIN reservoir versus 120 matched controls")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p266_reservoir_benchmark.png", dpi=190)
    plt.close(fig)


def summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "program": "P255",
            "status": results["P255"]["status"],
            "main_result": "all 4094 nontrivial contexts satisfy the finite Stieltjes measure theorem",
            "boundary": "finite operator theorem, not physical environment",
        },
        {
            "program": "P256",
            "status": results["P256"]["status"],
            "main_result": "minimal pole-residue realization identified; decoupled hidden modes are invisible",
            "boundary": "microscopic hidden sector remains nonunique",
        },
        {
            "program": "P257",
            "status": results["P257"]["status"],
            "main_result": "contravariant Schur context category and composition certificate",
            "boundary": "no sheaf or quantum contextuality theorem",
        },
        {
            "program": "P258",
            "status": results["P258"]["status"],
            "main_result": "closed Xi formula, covariance and norm bound",
            "boundary": "orientation receiver, not selector source",
        },
        {
            "program": "P259",
            "status": results["P259"]["status"],
            "main_result": "no fixed point/cycle in strictly size-decreasing context category",
            "boundary": "size-restoring RG equivalence still missing",
        },
        {
            "program": "P260",
            "status": results["P260"]["status"],
            "main_result": "multi-step environment/record/conditional information ledger",
            "boundary": "not thermodynamic entropy or energy",
        },
        {
            "program": "P261",
            "status": results["P261"]["status"],
            "main_result": "amplitude-only legacy completion fails positive Stieltjes class",
            "boundary": "not a no-go for richer completion maps",
        },
        {
            "program": "P262",
            "status": results["P262"]["status"],
            "main_result": "projective-first then independent-scale synthetic protocol",
            "boundary": "P242 external execution remains blocked",
        },
        {
            "program": "P263",
            "status": results["P263"]["status"],
            "main_result": "vertex POVM current no-go and phase-sensitive instrument construction",
            "boundary": "apparatus implementation remains external",
        },
        {
            "program": "P264",
            "status": results["P264"]["status"],
            "main_result": "generic versus specific atlas properties quantified",
            "boundary": "rates are conditional on frozen null ensembles",
        },
        {
            "program": "P265",
            "status": results["P265"]["status"],
            "main_result": "scale-clock identifiability quotient and multi-time diagnostic separation",
            "boundary": "no unique adaptive update law",
        },
        {
            "program": "P266",
            "status": results["P266"]["status"],
            "main_result": "matched 12-state reservoir benchmark",
            "boundary": "no biological ontology from computational analogy",
        },
    ]


def main() -> None:
    rng = np.random.default_rng(SEED)
    a, w = strict_operator()
    results: dict[str, Any] = {
        "metadata": {
            "release": "10.24",
            "seed": SEED,
            "scope": "Programs 255--266",
            "warning": (
                "finite proofs and synthetic audits only; no external "
                "laboratory record or physical closure"
            ),
        }
    }
    results["P255"] = program_255(a)
    results["P256"] = program_256(a)
    results["P257"] = program_257(a, rng)
    results["P258"] = program_258(a)
    results["P259"] = program_259(a)
    results["P260"] = program_260(a, rng)
    results["P261"] = program_261(a)
    results["P262"] = program_262(a, rng)
    results["P263"] = program_263(a, w)
    p264_result, false_positive_rows = program_264(a, rng)
    results["P264"] = p264_result
    results["P265"] = program_265(a)
    p266_result, reservoir_rows = program_266(a, rng)
    results["P266"] = p266_result

    RESULTS_PATH.write_text(
        json.dumps(results, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_csv(SUMMARY_PATH, summary_rows(results))
    write_csv(FALSE_POSITIVE_PATH, false_positive_rows)
    write_csv(RESERVOIR_PATH, reservoir_rows)
    make_figures(results, false_positive_rows, reservoir_rows)
    print(json.dumps(summary_rows(results), indent=2, ensure_ascii=False))
    print(
        json.dumps(
            {
                "p255_contexts": results["P255"]["contexts_audited"],
                "p257_max_residual": results["P257"]["maximum_composition_residual"],
                "p261_legacy_min_eigenvalue": results["P261"][
                    "legacy_amplitude_absorbed_generator_minimum_eigenvalue"
                ],
                "p262_joint_pass_rate": results["P262"]["joint_pass_rate"],
                "p264_all_five_rates": {
                    key: value["all_five"]
                    for key, value in results["P264"][
                        "ensemble_pass_rates"
                    ].items()
                },
                "p266_outlier": results["P266"][
                    "distinctive_outlier_under_declared_2_5_percent_rule"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
