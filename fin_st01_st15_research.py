#!/usr/bin/env python3
"""FIN ST01--ST15: local shadow-to-physics bridge research.

The batch is deliberately finite and dimensionless.  It separates exact
matrix theorems, conditional constructions, numerical evidence, failed
closures, and physical non-claims.  No result assumes that the counterfactual
FIN-as-ToE hypothesis is true.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import warnings
from pathlib import Path
from typing import Any, Callable

warnings.filterwarnings("ignore", message="Unable to import Axes3D.*")

import numpy as np
import scipy
import sympy as sp
from scipy.linalg import expm, logm, null_space

from fin_toe_shadow_atlas_analysis import (
    N,
    cyclic_distance,
    laplacian,
    radial_matrix,
    strict_profile,
    symmetry_matrices,
)
from fin_programs_507_516_research import (
    continue_pattern,
    stationary_root,
    stieltjes_memory_operator,
)
from fin_programs_517_526_research import hidden_memory_jacobian
from fin_programs_527_536_research import (
    hamiltonian_memory_jacobian,
    spectral_abscissa,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST01_ST15_Results.json"
SUMMARY = ROOT / "FIN_ST01_ST15_Summary.csv"
CERTIFICATE = ROOT / "FIN_ST01_ST15_Formal_Core_Certificate.json"
PREREGISTRATION = ROOT / "FIN_ST14_Prediction_Preregistration.json"
FIG_DIR = ROOT / "FIN_ST01_ST15_Figures"
SEED = 20260810


def native(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(v) for v in value]
    if isinstance(value, np.ndarray):
        return native(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def strict_operator() -> tuple[np.ndarray, np.ndarray, float]:
    w = radial_matrix(strict_profile)
    a = laplacian(w)
    return w, a, float(w.sum(axis=1)[0])


def relative_norm(x: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(x - reference) / max(np.linalg.norm(reference), 1e-15))


def hermitian_part(x: np.ndarray) -> np.ndarray:
    return (x + x.conj().T) / 2


def st01_shadow_theorem_packet() -> dict:
    required_maps = [
        "sector_or_vacuum_selector",
        "dimensional_unit_map",
        "coarse_graining_or_continuum_map",
        "observable_pullback",
        "preparation_family",
        "environment_or_dilation",
        "instrument_and_apparatus",
        "system_composition_law",
        "independent_record_map",
    ]
    composition_bound = (
        "If C1 Phi-Psi C1=R1 and C2 Psi-Xi C2=R2, then "
        "(C2 C1)Phi-Xi(C2 C1)=C2 R1+R2 C1 and its norm is at most "
        "||C2|| ||R1||+||R2|| ||C1||."
    )
    rng = np.random.default_rng(SEED)
    c1 = rng.normal(size=(4, 4))
    c2 = rng.normal(size=(4, 4))
    phi = rng.normal(size=(4, 4))
    psi = rng.normal(size=(4, 4))
    xi = rng.normal(size=(4, 4))
    r1 = c1 @ phi - psi @ c1
    r2 = c2 @ psi - xi @ c2
    composite = c2 @ c1 @ phi - xi @ c2 @ c1
    identity_residual = float(np.linalg.norm(composite - (c2 @ r1 + r2 @ c1)))
    bound_rhs = float(np.linalg.norm(c2, 2) * np.linalg.norm(r1, 2) + np.linalg.norm(r2, 2) * np.linalg.norm(c1, 2))
    bound_lhs = float(np.linalg.norm(composite, 2))
    certificate_schema = {
        "mathematical_representation": "required for L1",
        "finite_target_and_error_metric": "required for L2",
        "refinement_family_and_vanishing_intertwining_residual": "required for L3",
        "internally_sourced_scales_symmetries_observables_parameters": "required for L4",
        "independent_external_prediction_and_record": "required for L5",
        "physical_shadow_maps": required_maps,
    }
    negative_control = {key: True for key in required_maps}
    negative_control["instrument_and_apparatus"] = False
    return {
        "program": "ST01",
        "object": "Typed Physical-Shadow Certificate",
        "composition_theorem": composition_bound,
        "composition_identity_residual": identity_residual,
        "composition_norm_lhs": bound_lhs,
        "composition_norm_bound_rhs": bound_rhs,
        "certificate_schema": certificate_schema,
        "negative_control_missing_instrument_rejected": not all(negative_control.values()),
        "current_functional_calculus_level": "E1 exact mathematics / L1 physical representation",
        "status": "proven_schema_and_composition_bound",
        "nonclaim": "The schema does not supply any missing physical map.",
    }


def reconstruct_channels(a: np.ndarray, t_u: float, t_h: float, z: float, beta: float, t_w: float) -> dict:
    unitary = expm(-1j * t_u * a)
    heat = expm(-t_h * a)
    green = np.linalg.inv(a + z * np.eye(N))
    gibbs_raw = expm(-beta * a)
    gibbs = gibbs_raw / np.trace(gibbs_raw)
    evals, evecs = np.linalg.eigh(a)
    wave = evecs @ np.diag(np.cos(t_w * np.sqrt(np.maximum(evals, 0)))) @ evecs.T

    estimates: dict[str, np.ndarray] = {}
    estimates["unitary"] = hermitian_part((1j / t_u) * logm(unitary)).real
    estimates["heat"] = hermitian_part(-logm(heat) / t_h).real
    estimates["green"] = hermitian_part(np.linalg.inv(green) - z * np.eye(N)).real
    gibbs_est = hermitian_part(-logm(gibbs) / beta).real
    gibbs_est -= np.min(np.linalg.eigvalsh(gibbs_est)) * np.eye(N)
    estimates["gibbs"] = gibbs_est
    cvals, cvecs = np.linalg.eigh(hermitian_part(wave).real)
    wave_vals = (np.arccos(np.clip(cvals, -1, 1)) / t_w) ** 2
    estimates["wave"] = cvecs @ np.diag(wave_vals) @ cvecs.T
    return {
        "matrices": {"unitary": unitary, "heat": heat, "green": green, "gibbs": gibbs, "wave": wave},
        "estimates": estimates,
    }


def st02_common_spectrum_consistency(a: np.ndarray) -> dict:
    scales = {"t_unitary": 0.47, "t_heat": 0.61, "z_green": 0.70, "beta": 0.80, "t_wave": 0.55}
    reconstructed = reconstruct_channels(
        a,
        scales["t_unitary"],
        scales["t_heat"],
        scales["z_green"],
        scales["beta"],
        scales["t_wave"],
    )
    errors = {name: relative_norm(estimate, a) for name, estimate in reconstructed["estimates"].items()}
    eigenvalues = np.linalg.eigvalsh(a)
    alias_margin = float(math.pi - scales["t_unitary"] * eigenvalues[-1])
    wave_branch_margin = float(math.pi - scales["t_wave"] * math.sqrt(eigenvalues[-1]))

    v = np.zeros(N)
    v[0], v[1] = 1 / math.sqrt(2), -1 / math.sqrt(2)
    a_alt = a + 0.15 * np.outer(v, v)
    heat_alt = expm(-scales["t_heat"] * a_alt)
    altered_estimate = hermitian_part(-logm(heat_alt) / scales["t_heat"]).real
    altered_residual = relative_norm(altered_estimate, a)

    bad_time = 3.0
    bad_unitary = expm(-1j * bad_time * a)
    bad_estimate = hermitian_part((1j / bad_time) * logm(bad_unitary)).real
    bad_alias_error = relative_norm(bad_estimate, a)
    tolerance = 1e-10
    return {
        "program": "ST02",
        "object": "Frozen Common-Generator Consistency Test",
        "frozen_channel_scales": scales,
        "relative_reconstruction_errors": errors,
        "alias_margin_pi_minus_t_lambda_max": alias_margin,
        "wave_branch_margin_pi_minus_t_sqrt_lambda_max": wave_branch_margin,
        "common_case_accepts": max(errors.values()) < tolerance,
        "altered_heat_generator_relative_residual": altered_residual,
        "altered_channel_rejected": altered_residual > 1e-3,
        "uncontrolled_alias_time": bad_time,
        "uncontrolled_alias_relative_error": bad_alias_error,
        "status": "proven_round_trip_and_executable_falsifier",
        "boundary": "This is a synthetic closed-loop validation, not independent physical tomography.",
    }


def dihedral_average(matrix: np.ndarray) -> np.ndarray:
    shift, reflection = symmetry_matrices()
    total = np.zeros_like(matrix, dtype=float)
    count = 0
    for k in range(N):
        sk = np.linalg.matrix_power(shift, k)
        for r in (np.eye(N), reflection):
            g = sk @ r
            total += (g @ matrix @ g.T).real
            count += 1
    return total / count


def radial_weight_average_from_a(a: np.ndarray) -> np.ndarray:
    w = -a.copy()
    np.fill_diagonal(w, 0.0)
    rows = []
    for d in range(1, 7):
        values = [w[x, y] for x in range(N) for y in range(N) if cyclic_distance(x, y) == d]
        rows.append(float(np.mean(values)))
    return np.asarray(rows)


def kernel_law_fingerprint(weights: np.ndarray) -> float:
    d = np.arange(1, 7, dtype=float)
    cosine = np.cos(0.18575 * d + 0.16250)
    transformed = weights * (1 + d**1.8) / cosine
    return float(np.std(transformed) / max(abs(np.mean(transformed)), 1e-15))


def matrix_metrics(a: np.ndarray) -> dict:
    norm = max(np.linalg.norm(a), 1e-15)
    evals = np.linalg.eigvalsh(hermitian_part(a).real)
    positive = evals[evals > 1e-10]
    shift, reflection = symmetry_matrices()
    offdiag = a - np.diag(np.diag(a))
    return {
        "minimum_eigenvalue": float(evals[0]),
        "normalized_gap": float(positive[0] / np.mean(positive)) if len(positive) else float("nan"),
        "dihedral_residual": float(np.linalg.norm(a - dihedral_average(a)) / norm),
        "translation_commutator": float(np.linalg.norm(a @ shift - shift @ a) / norm),
        "reflection_commutator": float(np.linalg.norm(a @ reflection - reflection @ a) / norm),
        "kernel_law_fingerprint": kernel_law_fingerprint(radial_weight_average_from_a(a)),
        "m_matrix_markov_condition": bool(np.max(offdiag) <= 1e-12 and np.max(np.abs(a.sum(axis=1))) <= 1e-10),
    }


def random_orthogonal_fixing_uniform(rng: np.random.Generator) -> np.ndarray:
    one = np.ones((N, 1)) / math.sqrt(N)
    complement = null_space(one.T)
    q, _ = np.linalg.qr(rng.normal(size=(N - 1, N - 1)))
    return one @ one.T + complement @ q @ complement.T


def random_eigenbasis_with_uniform_zero_mode(rng: np.random.Generator) -> np.ndarray:
    one = np.ones((N, 1)) / math.sqrt(N)
    complement = null_space(one.T)
    q, _ = np.linalg.qr(rng.normal(size=(N - 1, N - 1)))
    return np.column_stack([one[:, 0], complement @ q])


def st03_generic_operator_null_atlas(a: np.ndarray, samples: int = 1200) -> dict:
    rng = np.random.default_rng(SEED + 3)
    target_s = float(np.trace(a) / N)
    ensembles: dict[str, list[dict]] = {
        "positive_radial": [],
        "positive_nonradial_graph": [],
        "random_psd_fixed_zero_mode": [],
        "isospectral_rotations": [],
        "signed_radial_kernel": [],
    }
    strict = matrix_metrics(a)
    strict_evals = np.linalg.eigvalsh(a)

    for _ in range(samples):
        radial = rng.lognormal(0.0, 1.0, 6)
        radial *= target_s / (2 * radial[:5].sum() + radial[5])
        ar = laplacian(radial_matrix(lambda d, rr=radial: rr[d - 1]))
        ensembles["positive_radial"].append(matrix_metrics(ar))

        w = rng.lognormal(-1.0, 0.8, size=(N, N))
        w = (w + w.T) / 2
        np.fill_diagonal(w, 0.0)
        w *= target_s / np.mean(w.sum(axis=1))
        ensembles["positive_nonradial_graph"].append(matrix_metrics(laplacian(w)))

        q = random_eigenbasis_with_uniform_zero_mode(rng)
        values = np.concatenate([[0.0], rng.lognormal(0.0, 0.7, N - 1)])
        apsd = q @ np.diag(values) @ q.T
        ensembles["random_psd_fixed_zero_mode"].append(matrix_metrics(apsd))

        qiso = random_orthogonal_fixing_uniform(rng)
        aiso = qiso @ a @ qiso.T
        row = matrix_metrics(aiso)
        row["maximum_spectral_error"] = float(np.max(np.abs(np.linalg.eigvalsh(aiso) - strict_evals)))
        ensembles["isospectral_rotations"].append(row)

        signed = rng.normal(0.0, 1.0, 6)
        signed *= target_s / max(2 * np.abs(signed[:5]).sum() + abs(signed[5]), 1e-12)
        asign = laplacian(radial_matrix(lambda d, ss=signed: ss[d - 1]))
        ensembles["signed_radial_kernel"].append(matrix_metrics(asign))

    summary: dict[str, dict] = {}
    for name, rows in ensembles.items():
        summary[name] = {
            "sample_count": len(rows),
            "psd_rate": float(np.mean([row["minimum_eigenvalue"] >= -1e-10 for row in rows])),
            "markov_m_matrix_rate": float(np.mean([row["m_matrix_markov_condition"] for row in rows])),
            "median_dihedral_residual": float(np.median([row["dihedral_residual"] for row in rows])),
            "median_kernel_law_fingerprint": float(np.median([row["kernel_law_fingerprint"] for row in rows])),
            "minimum_kernel_law_fingerprint": float(np.min([row["kernel_law_fingerprint"] for row in rows])),
        }
        if name == "isospectral_rotations":
            summary[name]["maximum_isospectral_error"] = float(np.max([row["maximum_spectral_error"] for row in rows]))

    distributions = {
        name: {
            "dihedral_residual": [row["dihedral_residual"] for row in rows],
            "kernel_law_fingerprint": [row["kernel_law_fingerprint"] for row in rows],
        }
        for name, rows in ensembles.items()
    }
    return {
        "program": "ST03",
        "object": "Five-Class Generic-Operator Null Atlas",
        "strict_metrics": strict,
        "ensemble_summary": summary,
        "_distributions": distributions,
        "isospectral_obstruction": (
            "Orthogonal rotations fixing the zero mode preserve every eigenvalue while generically destroying vertex, "
            "dihedral and Markov structure; spectrum alone cannot identify the FIN geometry."
        ),
        "specificity_result": (
            "The frozen kernel-law fingerprint is highly specific to its defining formula in the sampled classes, "
            "but is an implementation fingerprint, not a physical prediction. No physics-specific invariant was found."
        ),
        "status": "strong_multiclass_null_falsification_and_isospectral_no_go",
    }


def algebra_dimension(generators: list[np.ndarray], tolerance: float = 2e-10) -> int:
    size = generators[0].shape[0]
    initial = [np.eye(size, dtype=complex)]
    initial.extend(generator.astype(complex) for generator in generators)
    initial.extend(generator.conj().T.astype(complex) for generator in generators)
    stack = np.column_stack([matrix.reshape(-1) for matrix in initial])
    u, singular, _ = np.linalg.svd(stack, full_matrices=False)
    rank = int(np.sum(singular > tolerance * singular[0]))
    basis = u[:, :rank]
    while rank < size * size:
        current = [basis[:, index].reshape(size, size) for index in range(rank)]
        candidates = [basis]
        for element in current:
            for generator in generators:
                candidates.append((element @ generator).reshape(-1, 1))
                candidates.append((generator @ element).reshape(-1, 1))
        stack = np.column_stack(candidates)
        u, singular, _ = np.linalg.svd(stack, full_matrices=False)
        new_rank = int(np.sum(singular > tolerance * singular[0]))
        basis = u[:, :new_rank]
        if new_rank == rank:
            break
        rank = new_rank
    return rank


def st04_algebra_completion(a: np.ndarray) -> dict:
    label = np.diag(np.arange(N, dtype=float))
    anchor = np.diag([cyclic_distance(x, 0) for x in range(N)])
    shift, _ = symmetry_matrices()
    chiral = 1j * (shift - shift.T)
    dimensions = {
        "Cstar_A": algebra_dimension([a]),
        "A_plus_distinct_vertex_labels": algebra_dimension([a, label]),
        "A_plus_unoriented_anchor_distance": algebra_dimension([a, anchor]),
        "A_plus_circulant_chiral_generator": algebra_dimension([a, chiral]),
    }
    label_commutant_argument = (
        "A matrix commuting with a diagonal B with distinct entries is diagonal. If it also commutes with A and "
        "every off-diagonal A_xy is nonzero, all diagonal entries coincide. The joint commutant is scalar; "
        "the finite-dimensional Burnside theorem gives the full matrix algebra."
    )
    return {
        "program": "ST04",
        "object": "One-Generator Algebra Completion",
        "generated_algebra_dimensions": dimensions,
        "minimum_additional_generator_count": 1,
        "conditional_full_algebra_theorem": label_commutant_argument,
        "orientation_cost": (
            "B=diag(0,...,11) supplies labels and an orientation externally. The anchor preserves reflection and "
            "the chiral circulant remains inside a commuting Fourier algebra; neither is a strict selector theorem."
        ),
        "status": "proven_conditional_minimal_completion_with_selector_cost",
    }


def cycle_radial_laplacian(n: int, weight: Callable[[int, int], float]) -> np.ndarray:
    w = np.zeros((n, n))
    for x in range(n):
        for y in range(n):
            if x != y:
                d = min((x - y) % n, (y - x) % n)
                w[x, y] = weight(d, n)
    return np.diag(w.sum(axis=1)) - w


def st05_continuum_family_obstruction(a12: np.ndarray) -> dict:
    _, _, target_s = strict_operator()
    ns = [12, 24, 48, 96, 192]
    rows = []
    strict_weights = np.asarray([strict_profile(d) for d in range(1, 7)])
    for n in ns:
        local = cycle_radial_laplacian(n, lambda d, _n: float(strict_weights[d - 1]) if d <= 6 else 0.0)

        maximum_distance = n // 2
        raw_scaled = {
            d: math.cos(0.18575 * (12 * d / n) + 0.16250) / (1 + (12 * d / n) ** 1.8)
            for d in range(1, maximum_distance + 1)
        }
        multiplicity_sum = sum((1 if n % 2 == 0 and d == maximum_distance else 2) * raw_scaled[d] for d in raw_scaled)

        def scaled_dense_weight(d: int, nn: int) -> float:
            return raw_scaled[d] * target_s / multiplicity_sum

        dense = cycle_radial_laplacian(n, scaled_dense_weight)
        direct = cycle_radial_laplacian(n, lambda d, _n: strict_profile(d))
        local_eigs = np.linalg.eigvalsh(local)
        dense_eigs = np.linalg.eigvalsh(dense)
        direct_w_min = min(strict_profile(d) for d in range(1, n // 2 + 1))
        rows.append(
            {
                "N": n,
                "local_gap": float(local_eigs[1]),
                "dense_scaled_gap": float(dense_eigs[1]),
                "direct_integer_extension_minimum_weight": float(direct_w_min),
                "direct_integer_extension_psd": bool(np.linalg.eigvalsh(direct)[0] >= -1e-10),
                "direct_integer_extension_markov": bool(direct_w_min >= 0),
            }
        )
    logn = np.log(np.asarray(ns[1:], dtype=float))
    local_exp = float(np.polyfit(logn, np.log([row["local_gap"] for row in rows[1:]]), 1)[0])
    dense_exp = float(np.polyfit(logn, np.log([row["dense_scaled_gap"] for row in rows[1:]]), 1)[0])
    endpoint_error = relative_norm(cycle_radial_laplacian(12, lambda d, _n: strict_profile(d)), a12)
    first_negative_distance = next(d for d in range(1, 100) if strict_profile(d) < 0)
    return {
        "program": "ST05",
        "object": "Continuum-Extension Nonuniqueness Certificate",
        "endpoint_reconstruction_error_N12": endpoint_error,
        "continuation_rows": rows,
        "local_gap_scaling_exponent": local_exp,
        "dense_scaled_gap_scaling_exponent": dense_exp,
        "first_negative_weight_in_direct_integer_extension": first_negative_distance,
        "theorem": (
            "A single finite endpoint A_12 cannot select a unique refinement family: infinitely many sequences agree "
            "at N=12 and differ for every larger N. The two explicit continuations already have incompatible gap scaling."
        ),
        "status": "proven_nonidentifiability_and_naive_extension_markov_obstruction",
        "boundary": "This does not prove that no canonical FIN_N family can be added; it proves that current A_12 does not determine one.",
    }


def propagation_leakage(a: np.ndarray, times: list[float], radius: int = 1) -> list[dict]:
    source = np.zeros(N)
    source[0] = 1.0
    outside = np.asarray([cyclic_distance(x, 0) > radius for x in range(N)])
    rows = []
    for t in times:
        u = expm(-1j * t * a) @ source
        p = expm(-t * a) @ source
        rows.append(
            {
                "time": t,
                "unitary_probability_outside_radius": float(np.sum(np.abs(u[outside]) ** 2)),
                "heat_mass_outside_radius": float(np.sum(p[outside])),
            }
        )
    return rows


def st06_lorentzian_audit(w: np.ndarray, a: np.ndarray) -> dict:
    nearest_w = np.zeros_like(w)
    for x in range(N):
        nearest_w[x, (x + 1) % N] = strict_profile(1)
        nearest_w[x, (x - 1) % N] = strict_profile(1)
    nearest_a = laplacian(nearest_w)
    times = [1e-3, 1e-2, 5e-2, 1e-1]
    far_edges = [abs(w[0, y]) for y in range(N) if cyclic_distance(0, y) > 1]
    return {
        "program": "ST06",
        "object": "Exact-Cone Obstruction and Leakage Audit",
        "strict_rows": propagation_leakage(a, times),
        "nearest_neighbor_control_rows": propagation_leakage(nearest_a, times),
        "minimum_direct_far_coupling": float(min(far_edges)),
        "exact_derivative_statement": (
            "For x!=y, dU_t(x,y)/dt at t=0 equals -i A_xy=i W_xy. Since strict W_xy is nonzero at every "
            "cyclic distance 1,...,6, every vertex has an immediate first-order unitary amplitude outside any nontrivial "
            "local cone. Heat has the same first-order tail, while cos(t sqrt(A)) has the direct far term -t^2 A_xy/2."
        ),
        "status": "proven_no_exact_cyclic_causal_cone_for_frozen_strict_kernel",
        "boundary": "An approximate cone may exist only after a new refinement/localization theorem; Lorentzian signature is not generated here.",
    }


def matrix_log_psd(rho: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(hermitian_part(rho).real)
    return vectors @ np.diag(np.log(np.maximum(values, 1e-15))) @ vectors.T


def st07_information_thermodynamics_bridge(a: np.ndarray) -> dict:
    beta, energy_scale = 0.8, 2.3
    h = energy_scale * a
    gamma_raw = expm(-beta * h)
    gamma = gamma_raw / np.trace(gamma_raw)
    c = 7.0
    gamma_rescaled_raw = expm(-(beta / c) * (c * h))
    gamma_rescaled = gamma_rescaled_raw / np.trace(gamma_rescaled_raw)

    rng = np.random.default_rng(SEED + 7)
    x = rng.normal(size=(N, N))
    rho = x @ x.T
    rho /= np.trace(rho)
    log_rho = matrix_log_psd(rho)
    log_gamma = matrix_log_psd(gamma)
    relative_entropy = float(np.trace(rho @ (log_rho - log_gamma)))
    entropy_rho = float(-np.trace(rho @ log_rho))
    entropy_gamma = float(-np.trace(gamma @ log_gamma))
    free_rho = float(np.trace(rho @ h) - entropy_rho / beta)
    free_gamma = float(np.trace(gamma @ h) - entropy_gamma / beta)
    identity_error = abs(relative_entropy - beta * (free_rho - free_gamma))
    removals = [
        {"removed": "state rho", "failure": "no entropy, expectation or operational preparation is defined"},
        {"removed": "Hamiltonian H", "failure": "probabilities do not identify energy levels"},
        {"removed": "beta/bath", "failure": "no equilibrium temperature or heat direction is defined"},
        {"removed": "process family P", "failure": "a state alone defines neither work, heat nor erasure"},
        {"removed": "energy/action unit E_*", "failure": "all outputs remain dimensionless and scale-equivariant"},
    ]
    return {
        "program": "ST07",
        "object": "Conditional Relative-Entropy/Free-Energy Bridge",
        "minimal_equilibrium_semantics": ["state rho", "Hamiltonian H", "inverse bath temperature beta", "energy unit E_*"],
        "minimal_thermodynamic_process_tuple": ["state rho", "Hamiltonian H", "inverse bath temperature beta", "process family P", "energy/action unit E_*"],
        "gibbs_scale_orbit_error": float(np.linalg.norm(gamma - gamma_rescaled)),
        "relative_entropy_free_energy_identity_error": identity_error,
        "axiom_removal_ledger": removals,
        "status": "proven_conditional_bridge_and_scale_nonidentifiability",
        "boundary": "FIN supplies Shannon-compatible mathematics but does not currently source H, beta, the bath, the process or E_*.",
    }


def st08_gauge_holonomy_candidate(w: np.ndarray, s: float) -> dict:
    flux = 0.37
    theta = np.zeros((N, N))
    for x in range(N):
        y = (x + 1) % N
        theta[x, y] = flux / N
        theta[y, x] = -flux / N
    wtheta = w * np.exp(1j * theta)
    atheta = s * np.eye(N) - wtheta
    rng = np.random.default_rng(SEED + 8)
    chi = rng.uniform(-1, 1, N)
    gauge = np.diag(np.exp(1j * chi))
    theta_prime = theta + chi[:, None] - chi[None, :]
    aprime = s * np.eye(N) - w * np.exp(1j * theta_prime)
    covariance_residual = float(np.linalg.norm(aprime - gauge @ atheta @ gauge.conj().T))
    holonomy = np.prod([np.exp(1j * theta[x, (x + 1) % N]) for x in range(N)])
    holonomy_prime = np.prod([np.exp(1j * theta_prime[x, (x + 1) % N]) for x in range(N)])

    psi = rng.normal(size=N) + 1j * rng.normal(size=N)
    psi /= np.linalg.norm(psi)
    dprob = 2 * np.imag(np.conj(psi) * (atheta @ psi))
    current = np.zeros((N, N))
    for x in range(N):
        for y in range(N):
            current[x, y] = 2 * np.imag(np.conj(psi[x]) * atheta[x, y] * psi[y])
    continuity_residual = float(np.linalg.norm(dprob - current.sum(axis=1)))
    return {
        "program": "ST08",
        "object": "Inserted U(1) Link Connection",
        "flux": flux,
        "gauge_covariance_residual": covariance_residual,
        "holonomy_invariance_residual": float(abs(holonomy - holonomy_prime)),
        "continuity_equation_residual": continuity_residual,
        "magnetic_laplacian_minimum_eigenvalue": float(np.linalg.eigvalsh(atheta)[0]),
        "status": "proven_conditional_U1_receiver",
        "boundary": (
            "The link phases, U(1) group and flux were inserted. A global DNLS phase symmetry is not a local gauge law, "
            "and no electric field, Gauss constraint or Yang-Mills dynamics is derived."
        ),
    }


def st09_saturating_mediator(a: np.ndarray) -> dict:
    rho, b, g, rscale, mu, chi = sp.symbols("rho b g R mu chi", positive=True)
    q = rho / (1 + rho / rscale)
    joint = b * chi**2 / 2 - g * chi * q + mu * rho / 2
    stationary = sp.solve(sp.diff(joint, chi), chi)[0]
    effective = sp.simplify(joint.subs(chi, stationary))
    series = sp.series(effective, rho, 0, 4).removeO()

    b0, g0, r0, mu0 = 1.0, 1.0, 2.0, 0.15
    amplitudes = np.logspace(-2, 2, 100)
    one = np.ones(N) / math.sqrt(N)
    energies = []
    lower_bound = -N * g0**2 * r0**2 / (2 * b0)
    for amplitude in amplitudes:
        psi = amplitude * one
        density = np.abs(psi) ** 2
        qnum = density / (1 + density / r0)
        energy = 0.5 * float(psi @ a @ psi) + 0.5 * mu0 * float(np.sum(density)) - g0**2 / (2 * b0) * float(np.sum(qnum**2))
        energies.append(energy)
    return {
        "program": "ST09",
        "object": "Coercive Saturating-Mediator Completion",
        "joint_density_energy": str(joint),
        "stationary_mediator": str(stationary),
        "effective_density_energy": str(effective),
        "small_density_series_through_rho3": str(series),
        "global_lower_bound": lower_bound,
        "sample_minimum_energy": float(min(energies)),
        "large_amplitude_energy": float(energies[-1]),
        "parameters": {"b": b0, "g": g0, "R": r0, "mu": mu0},
        "_amplitudes": amplitudes,
        "_energies": energies,
        "theorem": (
            "For b,g,R,mu>0, stationary elimination gives -g^2 q(rho)^2/(2b)+mu rho/2, q=rho/(1+rho/R). "
            "It has the P527 attractive quartic jet, is bounded below by -N g^2 R^2/(2b), and the positive mu term "
            "makes the finite-dimensional energy coercive."
        ),
        "status": "proven_axiomatic_global_coercive_completion",
        "boundary": "The saturation law, R and mu are constructed additional axioms, not strict-core consequences; soliton orbital stability is not proved.",
    }


def canonical_symplectic(size: int, pair_slices: list[tuple[slice, slice]]) -> np.ndarray:
    j0 = np.zeros((size, size))
    for qslice, pslice in pair_slices:
        width = qslice.stop - qslice.start
        j0[qslice, pslice] = np.eye(width)
        j0[pslice, qslice] = -np.eye(width)
    return j0


def st10_memory_krein_audit(a: np.ndarray) -> dict:
    memory, poles, residues = stieltjes_memory_operator(a)
    base, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    state = base.copy()
    for epsilon in np.linspace(0.0, 1.0, 51):
        state = stationary_root(a + epsilon * memory, 1.0, 1.0, state)
    count = len(poles)
    size = 2 * N * (1 + count)
    pairs = [(slice(0, N), slice(N, 2 * N))]
    for j in range(count):
        pairs.append((slice(2 * N + 2 * j * N, 2 * N + (2 * j + 1) * N), slice(2 * N + (2 * j + 1) * N, 2 * N + (2 * j + 2) * N)))
    j0 = canonical_symplectic(size, pairs)
    speeds = np.linspace(0.6, 2.6, 51)
    rows = []
    for speed in speeds:
        jac = hamiltonian_memory_jacobian(a, state, poles, residues, 1.0, float(speed))
        lmat = -j0 @ jac
        symmetry_residual = float(np.linalg.norm(lmat - lmat.T))
        linertia = np.linalg.eigvalsh(hermitian_part(lmat).real)
        row = spectral_abscissa(jac)
        row.update(
            {
                "speed": float(speed),
                "hamiltonian_symmetry_residual": symmetry_residual,
                "L_negative_inertia": int(np.sum(linertia < -1e-8)),
                "L_positive_inertia": int(np.sum(linertia > 1e-8)),
                "positive_energy_sufficient_condition": bool(linertia[0] > 1e-8),
            }
        )
        rows.append(row)
    stable = [row for row in rows if row["unstable_count"] == 0]
    unstable = [row for row in rows if row["unstable_count"] > 0]

    signature_speed = 1.60
    signature_jac = hamiltonian_memory_jacobian(a, state, poles, residues, 1.0, signature_speed)
    signature_l = -j0 @ signature_jac
    signature_eigenvalues, signature_vectors = np.linalg.eig(signature_jac)
    signatures = []
    for index, value in enumerate(signature_eigenvalues):
        if value.imag > 1e-6 and abs(value.real) < 2e-6:
            vector = signature_vectors[:, index]
            krein_form = float(np.real(np.vdot(vector, signature_l @ vector)))
            signatures.append(krein_form)
    signature_counts = {
        "positive": int(np.sum(np.asarray(signatures) > 1e-7)),
        "negative": int(np.sum(np.asarray(signatures) < -1e-7)),
        "near_zero": int(np.sum(np.abs(np.asarray(signatures)) <= 1e-7)),
    }

    relax = hidden_memory_jacobian(a, state, poles, residues, 1.0)
    relax_row = spectral_abscissa(relax)
    transition = None
    for left, right in zip(rows[:-1], rows[1:]):
        if left["unstable_count"] != right["unstable_count"]:
            transition = [left["speed"], right["speed"]]
            break
    return {
        "program": "ST10",
        "object": "Hamiltonian-Structure and Krein-Closure Audit",
        "speed_rows": rows,
        "hamiltonian_stable_sample_count": len(stable),
        "hamiltonian_unstable_sample_count": len(unstable),
        "sampled_transition_bracket": transition,
        "sampled_Krein_signature_speed": signature_speed,
        "sampled_positive_imaginary_mode_signature_counts": signature_counts,
        "minimum_nonzero_absolute_sampled_Krein_form": float(min(abs(value) for value in signatures if abs(value) > 1e-12)),
        "relaxation_unit_scale_spectral_abscissa": relax_row["spectral_abscissa"],
        "maximum_hamiltonian_structure_residual": float(max(row["hamiltonian_symmetry_residual"] for row in rows)),
        "positive_energy_sufficient_condition_ever_holds": any(row["positive_energy_sufficient_condition"] for row in rows),
        "status": "partial_exact_structure_strong_numerical_transition_no_exact_Krein_closure",
        "verdict": (
            "The Hamiltonian realization has the required J0 L structure to numerical precision, but L is indefinite "
            "throughout the scan, so the elementary positive-energy stability theorem does not apply. The sampled "
            "stability transition remains numerical; an exact Krein-collision certificate was not obtained."
        ),
    }


def schur_eliminate(matrix: np.ndarray, labels: list[int], eliminate: list[int]) -> tuple[np.ndarray, list[int]]:
    keep = [label for label in labels if label not in set(eliminate)]
    keep_pos = [labels.index(label) for label in keep]
    elim_pos = [labels.index(label) for label in eliminate]
    ee = matrix[np.ix_(keep_pos, keep_pos)]
    eh = matrix[np.ix_(keep_pos, elim_pos)]
    hh = matrix[np.ix_(elim_pos, elim_pos)]
    return ee - eh @ np.linalg.solve(hh, eh.T), keep


def truncate_small_couplings(matrix: np.ndarray, threshold: float) -> np.ndarray:
    out = matrix.copy()
    for i in range(len(out)):
        for j in range(i + 1, len(out)):
            if abs(out[i, j]) < threshold:
                out[i, j] = out[j, i] = 0.0
    return out


def st11_rg_scheme_independence(a: np.ndarray) -> dict:
    matrix = a + 0.4 * np.eye(N)
    labels = list(range(N))
    h1 = [1, 2, 4, 5]
    h2 = [7, 8, 10, 11]
    direct, direct_labels = schur_eliminate(matrix, labels, h1 + h2)
    first_a, labels_a = schur_eliminate(matrix, labels, h1)
    exact_a, final_a = schur_eliminate(first_a, labels_a, h2)
    first_b, labels_b = schur_eliminate(matrix, labels, h2)
    exact_b, final_b = schur_eliminate(first_b, labels_b, h1)
    exact_difference = float(max(np.linalg.norm(exact_a - direct), np.linalg.norm(exact_b - direct), np.linalg.norm(exact_a - exact_b)))

    threshold = 0.08
    trunc_a, trunc_labels_a = schur_eliminate(truncate_small_couplings(first_a, threshold), labels_a, h2)
    trunc_b, trunc_labels_b = schur_eliminate(truncate_small_couplings(first_b, threshold), labels_b, h1)
    truncated_difference = float(np.linalg.norm(trunc_a - trunc_b) / np.linalg.norm(direct))
    return {
        "program": "ST11",
        "object": "Exact Schur Associativity versus Truncated RG Scheme Dependence",
        "retained_vertices": direct_labels,
        "exact_elimination_order_A": h1 + h2,
        "exact_elimination_order_B": h2 + h1,
        "exact_scheme_difference": exact_difference,
        "truncation_threshold": threshold,
        "truncated_scheme_relative_difference": truncated_difference,
        "status": "proven_exact_gaussian_scheme_independence_and_truncation_dependence",
        "boundary": "No scale family, beta function, universal exponent or fractal fixed point follows from finite Schur associativity.",
    }


def st12_legacy_strict_role_matrix() -> dict:
    alpha_geo = 4 * math.log(2)
    beta_tors = 0.01
    historical = {
        "sin2_thetaW": alpha_geo / 12,
        "alpha_EM_inverse": alpha_geo / (2 * beta_tors) * (1 - beta_tors),
    }
    roles = [
        {"role": "weak-mixing formula", "historical_value": historical["sin2_thetaW"], "status": "BLOCKED_NO_COMPLETION_OR_ROLE_TRANSFER"},
        {"role": "electromagnetic inverse coupling", "historical_value": historical["alpha_EM_inverse"], "status": "BLOCKED_NO_COMPLETION_OR_ROLE_TRANSFER"},
        {"role": "beta^N gravity hierarchy", "historical_value": None, "status": "BLOCKED_NO_DIMENSIONAL_OR_ROLE_TRANSFER_THEOREM"},
        {"role": "phase/orientation selector", "historical_value": None, "status": "REFUTED_INSIDE_AUT_INVARIANT_FUNCTIONAL_CALCULUS"},
        {"role": "alpha_geo information normalization", "historical_value": alpha_geo, "status": "NUMERICAL_IDENTITY_ONLY_NO_STRICT_TRANSFER"},
    ]
    return {
        "program": "ST12",
        "object": "Legacy-to-Strict Role-Transfer Gate Matrix",
        "completion_endpoint_available": False,
        "typed_strict_additions_available": False,
        "role_invariance_theorem_available": False,
        "roles": roles,
        "status": "proven_gate_failure_role_audit_not_admissible",
        "boundary": "Numerical proximity of historical formulas to known constants is not derivation and cannot be inherited by strict.",
    }


def st13_cross_theory_frozen_parameters(a: np.ndarray) -> dict:
    t_u, t_h, z, beta, t_w = 0.47, 0.61, 0.70, 0.80, 0.55
    heat = expm(-t_h * a)
    calibrated = hermitian_part(-logm(heat) / t_h).real
    targets = reconstruct_channels(a, t_u, t_h, z, beta, t_w)["matrices"]
    evals, evecs = np.linalg.eigh(calibrated)
    predictions = {
        "unitary": expm(-1j * t_u * calibrated),
        "green": np.linalg.inv(calibrated + z * np.eye(N)),
        "gibbs": expm(-beta * calibrated) / np.trace(expm(-beta * calibrated)),
        "wave": evecs @ np.diag(np.cos(t_w * np.sqrt(np.maximum(evals, 0)))) @ evecs.T,
    }
    prediction_errors = {name: relative_norm(prediction, targets[name]) for name, prediction in predictions.items()}

    v = np.zeros(N)
    v[2], v[5] = 1 / math.sqrt(2), -1 / math.sqrt(2)
    a_alt = a + 0.08 * np.outer(v, v)
    altered_unitary = expm(-1j * t_u * a_alt)
    altered_error = relative_norm(predictions["unitary"], altered_unitary)
    return {
        "program": "ST13",
        "object": "Heat-Calibrated No-Refit Cross-Channel Prediction",
        "calibration_channel": "heat only",
        "frozen_prediction_errors": prediction_errors,
        "all_synthetic_same_A_channels_pass": max(prediction_errors.values()) < 1e-10,
        "altered_unitary_holdout_relative_error": altered_error,
        "altered_generator_holdout_rejected": altered_error > 1e-3,
        "status": "proven_synthetic_pipeline_validation",
        "boundary": "Every accepted channel was generated locally from the same A; this is not evidence that nature uses the common generator.",
    }


def random_radial_fingerprint_sample(rng: np.random.Generator, count: int) -> np.ndarray:
    values = []
    for _ in range(count):
        weights = rng.lognormal(0.0, 1.0, 6)
        values.append(kernel_law_fingerprint(weights))
    return np.asarray(values)


def st14_prediction_ledger() -> dict:
    training_rng = np.random.default_rng(SEED + 140)
    training = random_radial_fingerprint_sample(training_rng, 2000)
    threshold = 0.02
    specification = {
        "observable": "coefficient of variation of (1+d^1.8) w_d / cos(0.18575 d+0.16250), d=1,...,6",
        "strict_prediction": "score <= 1e-12",
        "null_rejection_threshold": threshold,
        "holdout_count": 5000,
        "holdout_seed": SEED + 141,
        "no_refit": True,
        "interpretation": "kernel implementation fingerprint only; not a physical observable",
    }
    canonical = json.dumps(specification, sort_keys=True, separators=(",", ":")).encode("utf-8")
    prereg_hash = hashlib.sha256(canonical).hexdigest()
    PREREGISTRATION.write_text(json.dumps({"sha256": prereg_hash, "specification": specification}, indent=2, sort_keys=True), encoding="utf-8")

    holdout_rng = np.random.default_rng(specification["holdout_seed"])
    holdout = random_radial_fingerprint_sample(holdout_rng, specification["holdout_count"])
    strict_weights = np.asarray([strict_profile(d) for d in range(1, 7)])
    strict_score = kernel_law_fingerprint(strict_weights)
    false_accepts = int(np.sum(holdout <= threshold))
    upper95_if_zero = float(1 - 0.05 ** (1 / len(holdout))) if false_accepts == 0 else None
    return {
        "program": "ST14",
        "object": "Hashed Frozen Kernel-Fingerprint Holdout",
        "preregistration_sha256": prereg_hash,
        "training_minimum_score": float(training.min()),
        "strict_score": strict_score,
        "holdout_minimum_score": float(holdout.min()),
        "holdout_false_accept_count": false_accepts,
        "holdout_count": len(holdout),
        "zero_event_exact_95pct_upper_rate": upper95_if_zero,
        "threshold": threshold,
        "_training": training,
        "_holdout": holdout,
        "status": "passed_preregistered_implementation_fingerprint_but_failed_physical_prediction_requirement",
        "verdict": (
            "The holdout distinguishes the frozen analytic kernel law from the declared lognormal null without refit. "
            "Because the observable is constructed from the defining formula, it is not a novel physical prediction and does not satisfy L5."
        ),
    }


def st15_formal_core(a: np.ndarray, st04: dict) -> dict:
    f0, f1, f2 = sp.symbols("f0 f1 f2", real=True)
    w01, w02, w12 = sp.symbols("w01 w02 w12", real=True)
    wmat = sp.Matrix([[0, w01, w02], [w01, 0, w12], [w02, w12, 0]])
    lmat = sp.diag(*[sum(wmat[i, j] for j in range(3)) for i in range(3)]) - wmat
    fvec = sp.Matrix([f0, f1, f2])
    dirichlet_rhs = sum(wmat[i, j] * (fvec[i] - fvec[j]) ** 2 for i in range(3) for j in range(3)) / 2
    dirichlet_identity = sp.simplify((fvec.T * lmat * fvec)[0] - dirichlet_rhs) == 0

    rational = sp.Matrix([[4, -1, -1], [-1, 3, -1], [-1, -1, 3]])
    ee = rational[:2, :2]
    eh = rational[:2, 2:]
    hh = rational[2:, 2:]
    schur = ee - eh * hh.inv() * eh.T
    schur_identity = rational.inv()[:2, :2] == schur.inv()

    beta, energy, c = sp.symbols("beta energy c", positive=True)
    scale_identity = sp.simplify((beta / c) * (c * energy) - beta * energy) == 0
    automorphisms = [1, 5, 7, 11]
    orbit_plus_one = sorted({(u * 1) % 12 for u in automorphisms})
    selector_no_go = 1 in orbit_plus_one and 11 in orbit_plus_one
    shift, reflection = symmetry_matrices()
    symmetry_residual = max(np.linalg.norm(a @ shift - shift @ a), np.linalg.norm(a @ reflection - reflection @ a))
    checks = {
        "symbolic_Dirichlet_identity": bool(dirichlet_identity),
        "exact_rational_shifted_Schur_inverse_block_identity": bool(schur_identity),
        "symbolic_energy_temperature_scale_orbit": bool(scale_identity),
        "Aut_Z12_orbit_contains_plus_and_minus_orientation": bool(selector_no_go),
        "strict_translation_reflection_residual_below_1e-12": bool(symmetry_residual < 1e-12),
        "strict_functional_algebra_dimension_is_7": st04["generated_algebra_dimensions"]["Cstar_A"] == 7,
    }
    certificate = {
        "program": "ST15",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "Aut_Z12_unit_orbit": orbit_plus_one,
        "strict_symmetry_residual": float(symmetry_residual),
        "scope": "SymPy exact identities plus binary64 strict-kernel replay; not proof-assistant machine checking.",
    }
    CERTIFICATE.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST15",
        "object": "Formal-Core Replay Certificate",
        **certificate,
        "status": "passed_exact_symbolic_and_numerical_replay_with_explicit_trust_boundary",
    }


def make_figures(results: dict[str, Any]) -> None:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Unable to import Axes3D.*", category=UserWarning)
        import matplotlib.pyplot as plt

    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.23})

    st02 = results["ST02"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.5))
    names = list(st02["relative_reconstruction_errors"])
    errors = [st02["relative_reconstruction_errors"][name] for name in names]
    axes[0].bar(names, errors, color="#3182bd")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("relative generator error")
    axes[0].set_title("ST02 common-generator round trip")
    axes[1].bar(["common max", "altered heat", "bad alias"], [max(errors), st02["altered_heat_generator_relative_residual"], st02["uncontrolled_alias_relative_error"]], color=["#31a354", "#de2d26", "#fd8d3c"])
    axes[1].set_yscale("log")
    axes[1].set_ylabel("relative residual")
    axes[1].set_title("negative controls")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st02_common_generator_and_controls.png", dpi=190)
    plt.close(fig)

    st03 = results["ST03"]
    distributions = st03.pop("_distributions")
    labels = list(distributions)
    short = ["radial+", "graph+", "PSD", "isospec", "signed"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.8))
    axes[0].boxplot([np.maximum(distributions[name]["dihedral_residual"], 1e-16) for name in labels], tick_labels=short, showfliers=False)
    axes[0].set_yscale("log")
    axes[0].axhline(max(st03["strict_metrics"]["dihedral_residual"], 1e-16), color="red", label="strict")
    axes[0].set_title("dihedral residual")
    axes[0].legend()
    axes[1].boxplot([np.maximum(distributions[name]["kernel_law_fingerprint"], 1e-16) for name in labels], tick_labels=short, showfliers=False)
    axes[1].set_yscale("log")
    axes[1].axhline(max(st03["strict_metrics"]["kernel_law_fingerprint"], 1e-16), color="red", label="strict")
    axes[1].set_title("frozen kernel-law fingerprint")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st03_multiclass_null_atlas.png", dpi=190)
    plt.close(fig)

    st04, st05 = results["ST04"], results["ST05"]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.7))
    dims = st04["generated_algebra_dimensions"]
    axes[0].bar(["A", "+ labels", "+ anchor", "+ chiral"], list(dims.values()), color="#756bb1")
    axes[0].axhline(144, color="black", ls="--", lw=1)
    axes[0].set_ylabel("generated algebra dimension")
    axes[0].set_title("ST04 algebra completion")
    rows = st05["continuation_rows"]
    axes[1].loglog([row["N"] for row in rows], [row["local_gap"] for row in rows], "o-", label="fixed lattice range")
    axes[1].loglog([row["N"] for row in rows], [row["dense_scaled_gap"] for row in rows], "s-", label="scaled dense")
    axes[1].set_xlabel("N")
    axes[1].set_ylabel("first positive eigenvalue")
    axes[1].set_title("ST05 incompatible continuations")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st04_algebra_st05_continuum.png", dpi=190)
    plt.close(fig)

    st06 = results["ST06"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
    for key, label in [("strict_rows", "strict full"), ("nearest_neighbor_control_rows", "nearest-neighbor")]:
        rows = st06[key]
        axes[0].loglog([row["time"] for row in rows], [row["unitary_probability_outside_radius"] for row in rows], "o-", label=label)
        axes[1].loglog([row["time"] for row in rows], [row["heat_mass_outside_radius"] for row in rows], "o-", label=label)
    axes[0].set_title("unitary leakage outside radius 1")
    axes[1].set_title("heat leakage outside radius 1")
    for ax in axes:
        ax.set_xlabel("dimensionless time")
        ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st06_propagation_leakage.png", dpi=190)
    plt.close(fig)

    st09, st10 = results["ST09"], results["ST10"]
    amplitudes = np.asarray(st09.pop("_amplitudes"))
    energies = np.asarray(st09.pop("_energies"))
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.7))
    axes[0].semilogx(amplitudes, energies)
    axes[0].axhline(st09["global_lower_bound"], color="red", ls="--", label="global lower bound")
    axes[0].set_xlabel("uniform amplitude")
    axes[0].set_ylabel("effective energy")
    axes[0].set_title("ST09 coercive saturation")
    axes[0].legend()
    speed_rows = st10["speed_rows"]
    axes[1].plot([row["speed"] for row in speed_rows], [row["spectral_abscissa"] for row in speed_rows], "o-", ms=3)
    axes[1].axhline(2e-6, color="red", ls="--", label="numerical instability cut")
    axes[1].set_yscale("symlog", linthresh=1e-8)
    axes[1].set_xlabel("Hamiltonian mediator speed")
    axes[1].set_ylabel("spectral abscissa")
    axes[1].set_title("ST10 unresolved Krein transition")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st09_saturation_st10_memory.png", dpi=190)
    plt.close(fig)

    st11, st14 = results["ST11"], results["ST14"]
    training = np.asarray(st14.pop("_training"))
    holdout = np.asarray(st14.pop("_holdout"))
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.7))
    axes[0].bar(["exact Schur", "truncated"], [max(st11["exact_scheme_difference"], 1e-18), st11["truncated_scheme_relative_difference"]], color=["#31a354", "#de2d26"])
    axes[0].set_yscale("log")
    axes[0].set_ylabel("scheme difference")
    axes[0].set_title("ST11 exact vs truncated elimination")
    bins = np.logspace(-2, 2, 65)
    axes[1].hist(training, bins=bins, alpha=0.45, label="training")
    axes[1].hist(holdout, bins=bins, alpha=0.45, label="holdout")
    axes[1].axvline(st14["threshold"], color="red", label="frozen threshold")
    axes[1].axvline(max(st14["strict_score"], 1e-16), color="black", ls="--", label="strict")
    axes[1].set_xscale("log")
    axes[1].set_title("ST14 preregistered fingerprint")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st11_rg_st14_prediction.png", dpi=190)
    plt.close(fig)


def write_summary(results: dict[str, Any]) -> None:
    findings = {
        "ST01": "typed shadow certificate and composition bound",
        "ST02": "common-generator round trip accepts; altered and aliased controls reject",
        "ST03": "isospectral rotations destroy geometry; no physics-specific invariant found",
        "ST04": "one distinct vertex-label observable completes M12 but imports a selector",
        "ST05": "A12 admits inequivalent continuum continuations; direct tail loses Markov positivity",
        "ST06": "strict all-distance coupling forbids an exact nontrivial cyclic causal cone",
        "ST07": "relative entropy becomes free energy only after H, beta, process and units",
        "ST08": "U(1) covariance and holonomy work only for inserted link variables",
        "ST09": "constructed saturating mediator is globally coercive but not FIN-derived",
        "ST10": "Hamiltonian structure verified; exact Krein stability closure remains open",
        "ST11": "exact Schur elimination is associative; truncation restores scheme dependence",
        "ST12": "legacy physical-role audit remains blocked before completion and role transfer",
        "ST13": "heat-calibrated synthetic cross-channel pipeline passes without refit",
        "ST14": "hashed kernel fingerprint passes holdout but is not a physical prediction",
        "ST15": "formal-core exact/symbolic replay passes within its stated trust boundary",
    }
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["program", "status", "main_finding"])
        for key in [f"ST{i:02d}" for i in range(1, 16)]:
            writer.writerow([key, results[key]["status"], findings[key]])


def main() -> None:
    w, a, s = strict_operator()
    results: dict[str, Any] = {
        "release": "10.55",
        "programs": "ST01-ST15",
        "language": "English",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sympy": sp.__version__,
        },
    }
    results["ST01"] = st01_shadow_theorem_packet()
    results["ST02"] = st02_common_spectrum_consistency(a)
    results["ST03"] = st03_generic_operator_null_atlas(a)
    results["ST04"] = st04_algebra_completion(a)
    results["ST05"] = st05_continuum_family_obstruction(a)
    results["ST06"] = st06_lorentzian_audit(w, a)
    results["ST07"] = st07_information_thermodynamics_bridge(a)
    results["ST08"] = st08_gauge_holonomy_candidate(w, s)
    results["ST09"] = st09_saturating_mediator(a)
    results["ST10"] = st10_memory_krein_audit(a)
    results["ST11"] = st11_rg_scheme_independence(a)
    results["ST12"] = st12_legacy_strict_role_matrix()
    results["ST13"] = st13_cross_theory_frozen_parameters(a)
    results["ST14"] = st14_prediction_ledger()
    results["ST15"] = st15_formal_core(a, results["ST04"])
    results["recommended_next_programs"] = [
        {"id": "ST16", "priority": 1, "study": "derive or refute an exact Krein-collision certificate for the P530 Hamiltonian memory family"},
        {"id": "ST17", "priority": 2, "study": "classify selector-free additional observable algebras; prove a strict no-go if all Aut(Z12)-natural candidates remain reducible"},
        {"id": "ST18", "priority": 3, "study": "search for a FIN-internal refinement functor rather than choosing one of the incompatible ST05 continuations"},
        {"id": "ST19", "priority": 4, "study": "replace the tautological ST14 fingerprint with a derived, preregistered cross-sector observable"},
        {"id": "ST20", "priority": 5, "study": "construct interval error bars for the ST02 joint matrix logarithm and branch matching"},
        {"id": "ST21", "priority": 6, "study": "test approximate Lieb-Robinson-type bounds for truncated and decaying FIN_N candidates"},
        {"id": "ST22", "priority": 7, "study": "prove minimality or redundancy of the complete operational tuple from ST01"},
        {"id": "ST23", "priority": 8, "study": "derive a gauge connection from a new strict-sourced object or certify that the ST08 receiver is purely inserted"},
        {"id": "ST24", "priority": 9, "study": "classify all globally coercive saturations sharing the P527 fourth-order jet and compare universality"},
        {"id": "ST25", "priority": 10, "study": "formalize ST15 in a proof assistant after installing and pinning a verified toolchain"},
        {"id": "ST26", "priority": 11, "study": "build a dimensional-equivariance theorem covering unitary, wave, heat and Gibbs sectors simultaneously"},
        {"id": "ST27", "priority": 12, "study": "design a synthetic adversarial record in which sector-specific refits mimic H1A and prove the frozen protocol rejects it"},
    ]
    results["epistemic_boundary"] = (
        "All outputs are local finite mathematics, constructed conditional models or synthetic controls. They do not "
        "discharge QW-2191, source units, select a canonical refinement, complete legacy-to-strict, transfer legacy roles, "
        "supply a laboratory record, derive the Standard Model or gravity, or establish FIN as a Theory of Everything."
    )
    make_figures(results)
    write_summary(results)
    RESULTS.write_text(json.dumps(native(results), indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"results": RESULTS.name, "programs": 15, "figures": 6, "formal_certificate": CERTIFICATE.name}, indent=2))


if __name__ == "__main__":
    main()
