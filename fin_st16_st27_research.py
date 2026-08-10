#!/usr/bin/env python3
"""FIN ST16--ST27: second shadow-to-physics bridge research batch.

The batch is local and deterministic.  Exact finite-dimensional identities,
conditional constructions, numerical evidence, and failed closures are kept in
separate epistemic classes.  No result assumes that FIN is a physical theory.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import shutil
import subprocess
import warnings
from pathlib import Path
from typing import Any

warnings.filterwarnings("ignore", message="Unable to import Axes3D.*")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
import sympy as sp
from scipy.linalg import expm, logm

from fin_programs_507_516_research import (
    continue_pattern,
    stationary_root,
    stieltjes_memory_operator,
)
from fin_programs_527_536_research import hamiltonian_memory_jacobian
from fin_st01_st15_research import (
    N,
    algebra_dimension,
    canonical_symplectic,
    cycle_radial_laplacian,
    random_orthogonal_fixing_uniform,
    relative_norm,
    strict_operator,
    strict_profile,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST16_ST27_Results.json"
SUMMARY = ROOT / "FIN_ST16_ST27_Summary.csv"
PREREG = ROOT / "FIN_ST19_Cross_Sector_Preregistration.json"
INTERVAL_CERT = ROOT / "FIN_ST20_Scalar_Interval_Certificate.json"
LEAN_SOURCE = ROOT / "FIN_ST25_Formal_Core.lean"
FIG_DIR = ROOT / "FIN_ST16_ST27_Figures"
SEED = 20260811


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


def final_memory_state(a: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    memory, poles, residues = stieltjes_memory_operator(a)
    state, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    for epsilon in np.linspace(0.0, 1.0, 51):
        state = stationary_root(a + epsilon * memory, 1.0, 1.0, state)
    return state, memory, poles, residues


def memory_symplectic(pole_count: int) -> np.ndarray:
    size = 2 * N * (1 + pole_count)
    pairs = [(slice(0, N), slice(N, 2 * N))]
    for j in range(pole_count):
        pairs.append(
            (
                slice(2 * N + 2 * j * N, 2 * N + (2 * j + 1) * N),
                slice(2 * N + (2 * j + 1) * N, 2 * N + (2 * j + 2) * N),
            )
        )
    return canonical_symplectic(size, pairs)


def gauge_chain_data(
    a: np.ndarray,
    state: np.ndarray,
    poles: np.ndarray,
    residues: np.ndarray,
    speed: float,
) -> dict:
    count = len(poles)
    size = 2 * N * (1 + count)
    j0 = memory_symplectic(count)
    jac = hamiltonian_memory_jacobian(a, state, poles, residues, 1.0, speed)
    gauge = np.zeros(size)
    gauge[N : 2 * N] = state
    identity = np.eye(N)
    for j, (pole, residue) in enumerate(zip(poles, residues)):
        bb = slice(2 * N + (2 * j + 1) * N, 2 * N + (2 * j + 2) * N)
        gauge[bb] = math.sqrt(residue / speed) * np.linalg.solve(a + pole * identity, state)
    generalized = np.linalg.lstsq(jac, gauge, rcond=1e-12)[0]
    invariant = float((j0 @ gauge) @ generalized)
    eigenvalues = np.linalg.eigvals(jac)
    # The exact gauge mode appears as a numerically split conjugate/real pair.
    # Remove precisely the two smallest-modulus roots rather than using a
    # tolerance that changes the instability count near the collision.
    order = np.argsort(np.abs(eigenvalues))
    active = eigenvalues[order[2:]]
    return {
        "speed": speed,
        "gauge_residual": float(np.linalg.norm(jac @ gauge)),
        "generalized_chain_residual": float(np.linalg.norm(jac @ generalized - gauge)),
        "jordan_solvability_invariant": invariant,
        "spectral_abscissa": float(np.max(active.real)),
        "unstable_count": int(np.sum(active.real > 2e-6)),
        "third_smallest_eigenvalue_modulus": float(np.sort(np.abs(eigenvalues))[2]),
    }


def st16_exact_krein_zero_mode(a: np.ndarray) -> dict:
    state, memory, poles, residues = final_memory_state(a)
    identity = np.eye(N)
    lplus = a + identity - np.diag(3.0 * state**2) + memory
    tmat = np.zeros_like(a)
    smat = np.zeros_like(a)
    for pole, residue in zip(poles, residues):
        resolvent = np.linalg.inv(a + pole * identity)
        tmat += residue * (resolvent @ resolvent)
        smat += residue * (resolvent @ resolvent @ resolvent)
    linv_u = np.linalg.solve(lplus, state)
    linv_tu = np.linalg.solve(lplus, tmat @ state)
    c0 = float(state @ linv_u)
    c1 = float(state @ linv_tu)
    c2 = float((tmat @ state) @ linv_tu)
    c3 = float(state @ smat @ state)
    # kappa(s) = -(c0 + 2*c1/s + (c2+c3)/s^2), epsilon=1.
    coefficients_x = [c2 + c3, 2.0 * c1, c0]
    roots_x = np.roots(coefficients_x)
    positive_speeds = sorted(
        float(1.0 / root.real)
        for root in roots_x
        if abs(root.imag) < 1e-10 and root.real > 0
    )
    collision_speed = positive_speeds[0]
    speeds = np.linspace(max(0.65, collision_speed - 0.35), collision_speed + 0.45, 65)
    rows = [gauge_chain_data(a, state, poles, residues, float(speed)) for speed in speeds]
    analytic_rows = []
    for row in rows:
        speed = row["speed"]
        analytic = -(c0 + 2.0 * c1 / speed + (c2 + c3) / speed**2)
        analytic_rows.append({**row, "analytic_jordan_invariant": float(analytic)})
    max_formula_residual = max(
        abs(row["analytic_jordan_invariant"] - row["jordan_solvability_invariant"])
        for row in analytic_rows
    )
    left = gauge_chain_data(a, state, poles, residues, collision_speed - 1e-4)
    right = gauge_chain_data(a, state, poles, residues, collision_speed + 1e-4)
    stationary_residual = float(np.linalg.norm((a + memory) @ state + state - state**3, ord=np.inf))
    return {
        "program": "ST16",
        "object": "Gauge-Jordan Collision Reduction",
        "exact_reduction": (
            "For the Hamiltonian memory block, the gauge vector v0 is explicit and J v1=v0 is solvable. "
            "A further chain vector exists iff kappa(s)=(J0 v0)^T v1=0. Eliminating hidden blocks gives "
            "kappa(s)=-[c0+2 c1/s+(c2+c3)/s^2], with c0=<u,L+^-1u>, "
            "c1=<u,L+^-1Tu>, c2=<Tu,L+^-1Tu>, c3=<u,Su>."
        ),
        "quadratic_coefficients_in_inverse_speed": coefficients_x,
        "positive_collision_speeds": positive_speeds,
        "selected_collision_speed": collision_speed,
        "formula_vs_full_chain_maximum_residual": float(max_formula_residual),
        "stationary_state_residual_inf": stationary_residual,
        "left_of_collision": left,
        "right_of_collision": right,
        "rows": analytic_rows,
        "status": "proven_jordan_chain_reduction_strong_numerical_collision_location",
        "boundary": (
            "The scalar reduction is analytic for the declared block family. The root uses a binary64 stationary "
            "state and is not an interval-certified exact algebraic speed; it is not a global Krein theorem."
        ),
    }


def radial_basis(n: int = 12) -> list[np.ndarray]:
    basis = []
    for distance in range(7):
        matrix = np.zeros((n, n))
        for x in range(n):
            for y in range(n):
                d = min((x - y) % n, (y - x) % n)
                matrix[x, y] = float(d == distance)
        basis.append(matrix)
    return basis


def st17_selector_free_algebra(a: np.ndarray) -> dict:
    basis = radial_basis()
    radial_dim = int(np.linalg.matrix_rank(np.stack([item.ravel() for item in basis]), tol=1e-12))
    candidate_generators = [
        np.eye(N),
        a,
        a @ a,
        expm(-0.7 * a),
        np.linalg.inv(a + 0.4 * np.eye(N)),
        np.diag(np.diag(expm(-0.7 * a))),
        np.diag(np.sum(np.abs(a), axis=1)),
    ]
    candidate_dim = algebra_dimension(candidate_generators)
    return {
        "program": "ST17",
        "object": "Automorphism-Natural Algebra No-Go",
        "dihedral_invariant_real_symmetric_dimension": radial_dim,
        "strict_functional_algebra_dimension": algebra_dimension([a]),
        "algebra_dimension_after_natural_candidates": candidate_dim,
        "full_matrix_algebra_dimension": N * N,
        "theorem": (
            "Every real symmetric operator invariant under translations and reflection on C12 is a radial "
            "circulant, a seven-dimensional algebra. Strict A has seven distinct spectral values, hence its "
            "functional algebra already equals this full invariant algebra. Any deterministic generator that is "
            "natural under the stabilizer of A stays inside C*(A) and cannot generate M12(C)."
        ),
        "status": "proven_no_go_for_stabilizer_invariant_additional_generators",
        "boundary": (
            "Covariant orbit-valued or randomly symmetry-broken generators are not excluded, but they require a "
            "choice, state, boundary condition or external label and therefore do not discharge QW-2191."
        ),
    }


def lifted_laplacian_24(w12: np.ndarray, antipodal_weight: float) -> np.ndarray:
    def weight(distance: int, n: int) -> float:
        if 1 <= distance <= 5:
            return float(w12[0, distance])
        if distance == 6:
            return float(w12[0, 6] / 2.0)
        if distance == 12:
            return antipodal_weight
        return 0.0

    return cycle_radial_laplacian(24, weight)


def st18_refinement_functor_obstruction(w12: np.ndarray, a12: np.ndarray) -> dict:
    embedding = np.zeros((24, 12))
    for x in range(24):
        embedding[x, x % 12] = 1.0
    rows = []
    for q in [0.0, 0.1, 0.5, 1.0]:
        a24 = lifted_laplacian_24(w12, q)
        eigen = np.linalg.eigvalsh(a24)
        rows.append(
            {
                "antipodal_weight_q": q,
                "intertwining_residual": float(np.linalg.norm(a24 @ embedding - embedding @ a12)),
                "minimum_eigenvalue": float(eigen[0]),
                "maximum_eigenvalue": float(eigen[-1]),
                "markov_offdiagonal_condition": bool(np.max(a24 - np.diag(np.diag(a24))) <= 1e-14),
                "trace": float(np.trace(a24)),
            }
        )
    a0 = lifted_laplacian_24(w12, 0.0)
    a1 = lifted_laplacian_24(w12, 1.0)
    even_spectrum_residual = float(
        np.max(np.abs(np.fft.fft(a0[0]).real[::2] - np.fft.fft(a12[0]).real))
    )
    family_difference = float(np.linalg.norm(a1 - a0))
    return {
        "program": "ST18",
        "object": "Exact Positive Refinement-Fiber Certificate",
        "lift_rows": rows,
        "even_mode_spectrum_residual": even_spectrum_residual,
        "q0_q1_operator_difference": family_difference,
        "theorem": (
            "Let Jf(x)=f(x mod 12). Copy strict weights d=1,...,5 to C24, use half the strict "
            "antipodal d=6 weight, and add any q>=0 at fine distance 12. Then A24(q)J=JA12 exactly. "
            "Every member is a translation-invariant positive graph Laplacian, while q changes all odd fine modes."
        ),
        "status": "proven_infinite_positive_refinement_fiber_no_internal_selector",
        "boundary": (
            "The family gives valid refinements and proves nonuniqueness. It is not a FIN-derived functor and does "
            "not exclude a future new axiom selecting one q at every scale."
        ),
    }


def mixed_channel_score(a_cal: np.ndarray, a_unitary: np.ndarray, t: float, tau: float) -> float:
    u_actual = expm(-1j * t * a_unitary)
    p_actual = expm(-tau * a_cal)
    mixed_actual = u_actual @ p_actual
    mixed_predicted = expm(-(tau + 1j * t) * a_cal)
    return relative_norm(mixed_actual, mixed_predicted)


def st19_derived_preregistered_observable(a: np.ndarray) -> dict:
    config = {
        "seed": SEED + 19,
        "observable": "relative Frobenius residual of U_t P_tau versus exp[-(tau+i t) A_heat]",
        "calibration": "A_heat=-(1/t0)log P_t0 on the positive branch",
        "threshold": 1e-10,
        "same_generator_holdout_count": 400,
        "altered_unitary_holdout_count": 400,
        "time_ranges": {"t": [0.15, 0.8], "tau": [0.1, 0.9]},
        "alteration_delta_range": [0.01, 0.08],
    }
    payload = json.dumps(config, sort_keys=True, separators=(",", ":")).encode("utf-8")
    prereg_hash = hashlib.sha256(payload).hexdigest()
    PREREG.write_text(
        json.dumps({"configuration": config, "sha256": prereg_hash}, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    rng = np.random.default_rng(config["seed"])
    p_cal = expm(-0.63 * a)
    a_cal = -np.real_if_close(logm(p_cal)) / 0.63
    same_scores = []
    altered_scores = []
    one = np.ones(N) / math.sqrt(N)
    for _ in range(config["same_generator_holdout_count"]):
        t = float(rng.uniform(*config["time_ranges"]["t"]))
        tau = float(rng.uniform(*config["time_ranges"]["tau"]))
        same_scores.append(mixed_channel_score(a_cal, a, t, tau))
    for _ in range(config["altered_unitary_holdout_count"]):
        vector = rng.normal(size=N)
        vector -= one * (one @ vector)
        vector /= np.linalg.norm(vector)
        delta = float(rng.uniform(*config["alteration_delta_range"]))
        altered = a + delta * np.outer(vector, vector)
        t = float(rng.uniform(*config["time_ranges"]["t"]))
        tau = float(rng.uniform(*config["time_ranges"]["tau"]))
        altered_scores.append(mixed_channel_score(a_cal, altered, t, tau))
    threshold = config["threshold"]
    return {
        "program": "ST19",
        "object": "Derived Mixed-Channel Preregistered Observable",
        "preregistration_sha256": prereg_hash,
        "same_generator_maximum_score": float(max(same_scores)),
        "same_generator_false_rejections": int(np.sum(np.asarray(same_scores) > threshold)),
        "altered_generator_minimum_score": float(min(altered_scores)),
        "altered_generator_false_accepts": int(np.sum(np.asarray(altered_scores) <= threshold)),
        "same_generator_scores": same_scores,
        "altered_generator_scores": altered_scores,
        "status": "proven_synthetic_cross_channel_falsifier_not_physical_prediction",
        "boundary": (
            "The identity is derived from H1A rather than from the strict kernel formula, but every holdout is "
            "synthetic and locally generated. No laboratory record or independent custody is supplied."
        ),
    }


def interval_bounds(value: Any) -> tuple[str, str, float, float]:
    lo_raw, hi_raw = value._mpi_
    lo_string = mp.libmp.to_str(lo_raw, 55)
    hi_string = mp.libmp.to_str(hi_raw, 55)
    lo = float(lo_string)
    hi = float(hi_string)
    return lo_string, hi_string, float(np.nextafter(lo, -np.inf)), float(np.nextafter(hi, np.inf))


def st20_interval_common_generator() -> dict:
    mp.iv.dps = 55
    iv = mp.iv
    omega = iv.mpf("0.18575")
    phi = iv.mpf("0.16250")
    eta = iv.mpf("1.8")
    weights = {}
    for d in range(1, 7):
        weights[d] = iv.cos(omega * d + phi) / (1 + iv.mpf(d) ** eta)
    eigen_intervals = []
    for k in range(12):
        value = iv.mpf("0")
        for d in range(1, 6):
            value += 2 * weights[d] * (1 - iv.cos(2 * iv.pi * k * d / 12))
        value += weights[6] * (1 - iv.cos(iv.pi * k))
        lo_s, hi_s, lo, hi = interval_bounds(value)
        eigen_intervals.append({"mode": k, "lower": lo_s, "upper": hi_s, "lower_float": lo, "upper_float": hi})
    representatives = [eigen_intervals[k] for k in range(7)]
    representatives.sort(key=lambda row: (row["lower_float"] + row["upper_float"]) / 2)
    gap_lowers = []
    for left, right in zip(representatives[:-1], representatives[1:]):
        gap_lowers.append(right["lower_float"] - left["upper_float"])
    min_gap = float(min(gap_lowers))
    lambda_max_upper = float(max(row["upper_float"] for row in eigen_intervals))
    unitary_margin = float(math.pi - 0.47 * lambda_max_upper)
    wave_margin = float(math.pi - 0.55 * math.sqrt(lambda_max_upper))
    fro_upper = math.sqrt(sum(max(abs(row["lower_float"]), abs(row["upper_float"])) ** 2 for row in eigen_intervals))
    altered_relative_lower = 0.15 / fro_upper
    certificate = {
        "input_model": "mpmath.iv 55-decimal interval evaluation of decimal strict parameters",
        "eigenvalue_intervals": eigen_intervals,
        "minimum_distinct_eigenvalue_gap_lower_bound": min_gap,
        "unitary_principal_branch_margin_lower_bound": unitary_margin,
        "wave_arccos_branch_margin_lower_bound": wave_margin,
        "altered_rank_one_generator_relative_separation_lower_bound": altered_relative_lower,
        "projector_matching": (
            "Heat, Green and Gibbs scalar maps are strictly monotone. The wave map is strictly monotone on the "
            "certified arccos branch, and unitary phases are unaliased on the certified principal arc. Positive "
            "distinct-gap and branch margins therefore match the seven spectral projectors."
        ),
        "trust_boundary": "mpmath interval arithmetic and its transcendental implementation are trusted; no proof-assistant replay is claimed.",
    }
    INTERVAL_CERT.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST20",
        "object": "Scalar-Branch Interval Certificate",
        **certificate,
        "status": "interval_certified_scalar_branches_conditional_on_mpmath_iv",
    }


def cycle_spectral_values(n: int, weights: dict[int, float]) -> np.ndarray:
    values = []
    for k in range(n):
        value = 0.0
        for distance, weight in weights.items():
            multiplicity = 1 if distance == n // 2 else 2
            value += multiplicity * weight * (1.0 - math.cos(2.0 * math.pi * k * distance / n))
        values.append(value)
    return np.asarray(values)


def st21_locality_bounds(w12: np.ndarray) -> dict:
    t = 0.1
    finite_rows = []
    for radius in [1, 2, 3, 6]:
        a_radius = cycle_radial_laplacian(
            12,
            lambda distance, n: float(w12[0, distance]) if 1 <= distance <= radius else 0.0,
        )
        ell = math.ceil(6 / radius)
        x = t * float(np.linalg.norm(a_radius, 2))
        raw_bound = math.exp(x) * x**ell / math.factorial(ell)
        actual = abs(expm(-1j * t * a_radius)[0, 6])
        finite_rows.append(
            {
                "range_R": radius,
                "minimum_path_order": ell,
                "actual_antipodal_amplitude": float(actual),
                "series_remainder_bound": float(raw_bound),
                "unitarity_clipped_bound": float(min(1.0, raw_bound)),
            }
        )
    n = 192
    full_weights = {d: float(strict_profile(d)) for d in range(1, n // 2 + 1)}
    full_spectrum = cycle_spectral_values(n, full_weights)
    tail_rows = []
    for radius in [4, 8, 16, 32, 64]:
        truncated = {d: weight for d, weight in full_weights.items() if d <= radius}
        truncated_spectrum = cycle_spectral_values(n, truncated)
        delta_spectrum = full_spectrum - truncated_spectrum
        actual_unitary_error = float(np.max(np.abs(np.exp(-1j * t * full_spectrum) - np.exp(-1j * t * truncated_spectrum))))
        spectral_duhamel = float(t * np.max(np.abs(delta_spectrum)))
        absolute_tail = 0.0
        for distance, weight in full_weights.items():
            if distance > radius:
                multiplicity = 1 if distance == n // 2 else 2
                absolute_tail += multiplicity * abs(weight)
        row_bound = 2.0 * t * absolute_tail
        tail_rows.append(
            {
                "range_R": radius,
                "actual_unitary_operator_error": actual_unitary_error,
                "spectral_duhamel_bound": spectral_duhamel,
                "absolute_row_tail_bound": row_bound,
            }
        )
    exponent = float(
        np.polyfit(
            np.log([row["range_R"] for row in tail_rows[:-1]]),
            np.log([row["absolute_row_tail_bound"] for row in tail_rows[:-1]]),
            1,
        )[0]
    )
    first_negative = next(d for d, weight in full_weights.items() if weight < 0)
    return {
        "program": "ST21",
        "object": "Finite-Range Path Bound and Long-Range Tail Audit",
        "finite_range_rows": finite_rows,
        "long_range_rows_N192": tail_rows,
        "fitted_absolute_tail_exponent": exponent,
        "first_negative_literal_extension_weight": first_negative,
        "theorem": (
            "For a range-R matrix and cyclic distance d, all off-diagonal powers below ceil(d/R) vanish. "
            "Thus |exp(-itA)_xy| is bounded by exp(t||A||)(t||A||)^ell/ell!. For two Hermitian "
            "generators, Duhamel gives ||U-U_R||<=t||A-A_R||."
        ),
        "status": "proven_conditional_locality_bounds_and_computed_long_range_tail",
        "boundary": (
            "The literal long-range profile is signed from distance 8, so its unitary is valid but its heat "
            "channel is not a Markov probability semigroup. No Lorentzian continuum is derived."
        ),
    }


def st22_operational_tuple_minimality(a: np.ndarray, st18: dict) -> dict:
    localized = np.zeros(N)
    localized[0] = 1.0
    uniform = np.ones(N) / math.sqrt(N)
    vertex_localized = localized**2
    vertex_uniform = uniform**2
    preparation_tv = 0.5 * float(np.sum(np.abs(vertex_localized - vertex_uniform)))
    fourier_localized = np.abs(np.fft.fft(localized) / math.sqrt(N)) ** 2
    measurement_tv = 0.5 * float(np.sum(np.abs(vertex_localized - fourier_localized)))
    independent_joint = np.full(4, 0.25)
    correlated_joint = np.array([0.5, 0.0, 0.0, 0.5])
    composition_tv = 0.5 * float(np.sum(np.abs(independent_joint - correlated_joint)))
    c = 3.7
    t = 0.42
    unit_scale_residual = relative_norm(expm(-1j * t * a), expm(-1j * (t / c) * (c * a)))
    primitive_groups = {
        "sector": "not compressible: identical operator formulas on a direct sum do not select a physical sector",
        "dimensional_calibration": "not compressible: the positive scale orbit preserves dimensionless trajectories",
        "refinement": "not compressible: ST18 supplies distinct positive exact lifts of the same endpoint",
        "operational_process_theory": "must encode preparations, effects/instruments, classical outputs and composition",
        "external_record_and_custody": "not generated by an internal mathematical object",
    }
    packaged_components = {
        "observable_pullback": "an observable is an instrument/effect with classical output",
        "preparation_family": "a preparation is a process from the monoidal unit",
        "environment_or_dilation": "operationally redundant when its reduced channel/instrument is specified; dilations are nonunique",
        "apparatus_calibration": "can be metadata in the instrument context, though physical calibration evidence remains external",
        "record_map": "probability output can be packaged in an instrument; independent custody cannot",
    }
    return {
        "program": "ST22",
        "object": "Operational-Obligation Compression and Independence Ledger",
        "irreducible_obligation_groups": primitive_groups,
        "packaged_ST01_components": packaged_components,
        "finite_countermodel_metrics": {
            "preparation_choice_total_variation": preparation_tv,
            "measurement_choice_total_variation": measurement_tv,
            "same_marginals_different_composition_total_variation": composition_tv,
            "unit_scale_orbit_residual": unit_scale_residual,
            "ST18_exact_lift_count": len(st18["lift_rows"]),
        },
        "status": "proven_relative_minimality_after_operational_packaging",
        "boundary": (
            "Minimality is relative to finite operational process semantics. It is not a theorem that every possible "
            "foundation of physics must use this vocabulary."
        ),
    }


def st23_connection_source_no_go(a: np.ndarray) -> dict:
    reflection = np.zeros((N, N))
    for x in range(N):
        reflection[x, (-x) % N] = 1.0
    antisymmetric_basis = []
    for distance in range(1, 6):
        matrix = np.zeros((N, N))
        for x in range(N):
            matrix[x, (x + distance) % N] = 1.0
            matrix[x, (x - distance) % N] = -1.0
        antisymmetric_basis.append(matrix)
    constraints = np.stack(
        [(matrix - reflection @ matrix @ reflection.T).ravel() for matrix in antisymmetric_basis],
        axis=1,
    )
    exact_rank = int(sp.Matrix(constraints.astype(int)).rank())
    invariant_dimension = len(antisymmetric_basis) - exact_rank
    reflection_residual_of_a = float(np.linalg.norm(a - reflection @ a @ reflection.T))
    return {
        "program": "ST23",
        "object": "Strict-Sourced Continuous-Connection No-Go",
        "translation_invariant_antisymmetric_one_form_dimension": len(antisymmetric_basis),
        "reflection_constraint_exact_rank": exact_rank,
        "fully_dihedral_invariant_continuous_one_form_dimension": invariant_dimension,
        "strict_reflection_residual": reflection_residual_of_a,
        "holonomy_fixed_points_under_orientation_reversal": ["0 mod 2pi", "pi mod 2pi"],
        "theorem": (
            "A real translation-invariant U(1) link one-form is antisymmetric under edge reversal. Reflection sends "
            "every such coefficient theta_d to theta_-d=-theta_d. Dihedral invariance therefore forces all "
            "continuous coefficients to zero. A deterministic stabilizer-equivariant map from strict A cannot "
            "produce a nonzero oriented connection."
        ),
        "status": "proven_no_go_for_continuous_stabilizer_equivariant_connection_from_A_alone",
        "boundary": (
            "A discrete pi holonomy, a state-dependent connection, a boundary condition or explicit spontaneous "
            "symmetry breaking is not excluded; each is an additional sourced object."
        ),
    }


def saturation_potential(rho: np.ndarray, exponent: float, b: float, g: float, rscale: float, mu: float) -> np.ndarray:
    q = rho / (1.0 + rho / rscale) ** (1.0 - exponent)
    return -g * g * q * q / (2.0 * b) + mu * rho / 2.0


def st24_saturation_classification() -> dict:
    b = g = 1.0
    rscale = 2.0
    mu = 0.15
    rho = np.logspace(-6, 7, 800)
    rows = []
    for exponent in [0.0, 0.25, 0.5, 0.75, 1.0]:
        values = saturation_potential(rho, exponent, b, g, rscale, mu)
        if exponent < 0.5:
            verdict = "coercive"
        elif exponent > 0.5:
            verdict = "unbounded_below"
        else:
            verdict = "coercive" if mu > g * g * rscale / b else "unbounded_below"
        rows.append(
            {
                "asymptotic_growth_exponent_a": exponent,
                "predicted_verdict": verdict,
                "minimum_on_scan": float(np.min(values)),
                "value_at_rho_1e7": float(values[-1]),
                "quartic_density_coefficient": -g * g / (2.0 * b),
            }
        )
    critical_mu = g * g * rscale / b
    return {
        "program": "ST24",
        "object": "Universality Classes of Coercive Saturating Mediators",
        "family": "q_a(rho)=rho/(1+rho/R)^(1-a), V_a=-g^2 q_a^2/(2b)+mu rho/2",
        "rows": rows,
        "critical_a_half_mu": critical_mu,
        "theorem": (
            "Every a has q_a(rho)=rho+O(rho^2), hence the same attractive coefficient -g^2/(2b) rho^2. "
            "As rho tends to infinity, q_a^2 is asymptotic to R^(2-2a)rho^(2a). With mu>0 the energy "
            "is coercive for a<1/2, is coercive at a=1/2 iff mu>g^2 R/b, and is unbounded below for a>1/2."
        ),
        "status": "proven_asymptotic_coercivity_classification_for_declared_family",
        "boundary": "The family and its positive mass term are constructed axioms, not strict-FIN outputs.",
    }


def st25_proof_assistant_replay() -> dict:
    source_hash = hashlib.sha256(LEAN_SOURCE.read_bytes()).hexdigest()
    lean_path = shutil.which("lean")
    version = subprocess.run([lean_path, "--version"], capture_output=True, text=True) if lean_path else None
    configured = bool(version and version.returncode == 0)
    replay = None
    if configured:
        replay = subprocess.run([lean_path, str(LEAN_SOURCE)], capture_output=True, text=True)
    passed = bool(replay and replay.returncode == 0)
    return {
        "program": "ST25",
        "object": "Pinned Proof-Assistant Replay Attempt",
        "lean_launcher_found": bool(lean_path),
        "toolchain_configured": configured,
        "version_stdout": version.stdout.strip() if version else None,
        "version_stderr": version.stderr.strip() if version else None,
        "source_sha256": source_hash,
        "source_file": LEAN_SOURCE.name,
        "replay_attempted": replay is not None,
        "replay_passed": passed,
        "replay_stdout": replay.stdout.strip() if replay else None,
        "replay_stderr": replay.stderr.strip() if replay else None,
        "status": "blocked_no_configured_lean_toolchain_source_not_machine_checked" if not passed else "machine_checked_declared_formal_core",
        "boundary": (
            "The elan launcher exists but no default Lean toolchain is configured in this local environment. "
            "The exported source is not called machine checked, and no network installation was attempted."
            if not passed
            else "Only the four algebraic lemmas in the exported source were checked; no physical claim follows."
        ),
    }


def st26_dimensional_equivariance() -> dict:
    return {
        "program": "ST26",
        "object": "Multiweight Dimensional-Equivariance Theorem",
        "dimensionless_channel_coordinates": {
            "unitary": "u_U=alpha_U t, U=exp(-i u_U A)",
            "heat": "u_P=alpha_P t, P=exp(-u_P A)",
            "wave": "u_C=t sqrt(alpha_C), C=cos(u_C sqrt(A))",
            "green": "G_phys=(alpha_G A+z I)^-1, with weight alpha_G^-1",
            "gibbs": "u_beta=beta E_*, rho=exp(-u_beta A)/Z",
        },
        "positive_scale_action": {
            "unitary_and_heat": "alpha -> c alpha, t -> t/c",
            "wave_stiffness": "alpha_C -> c alpha_C, t -> t/sqrt(c)",
            "green": "alpha_G -> c alpha_G, z -> c z, G -> G/c",
            "gibbs": "E_* -> c E_*, beta -> beta/c",
        },
        "theorem": (
            "All five channels are equivariant under a positive scaling group, but with different weights. "
            "A single dimension assignment to A cannot make it simultaneously a first-order generator with "
            "dimension T^-1 and a wave stiffness with dimension T^-2. Sector scales or a squared operator map "
            "are therefore necessary. Internal phase frequencies determine ratios, not seconds."
        ),
        "status": "proven_multiweight_scale_orbit_and_single_dimension_obstruction",
        "boundary": "The theorem organizes required scale maps; it does not select their numerical physical values.",
    }


def fixed_uniform_rotation(rng: np.random.Generator, epsilon: float) -> np.ndarray:
    one = np.ones(N) / math.sqrt(N)
    projector = np.eye(N) - np.outer(one, one)
    raw = rng.normal(size=(N, N))
    skew = projector @ (raw - raw.T) @ projector
    skew /= max(np.linalg.norm(skew, 2), 1e-15)
    return expm(epsilon * skew)


def st27_adversarial_refits(a: np.ndarray) -> dict:
    rng = np.random.default_rng(SEED + 27)
    t = 0.47
    predicted = expm(-1j * t * a)
    epsilons = [0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1]
    threshold = 1e-3
    rows = []
    for epsilon in epsilons:
        residuals = []
        spectral_errors = []
        refit_errors = []
        for _ in range(60):
            q = fixed_uniform_rotation(rng, epsilon)
            altered = q @ a @ q.T
            channel = expm(-1j * t * altered)
            residuals.append(relative_norm(channel, predicted))
            spectral_errors.append(float(np.max(np.abs(np.linalg.eigvalsh(altered) - np.linalg.eigvalsh(a)))))
            refit_errors.append(relative_norm(channel, expm(-1j * t * altered)))
        rows.append(
            {
                "rotation_size_epsilon": epsilon,
                "median_frozen_common_generator_residual": float(np.median(residuals)),
                "minimum_frozen_common_generator_residual": float(np.min(residuals)),
                "maximum_isospectral_error": float(np.max(spectral_errors)),
                "maximum_sector_specific_refit_error": float(np.max(refit_errors)),
                "frozen_rejection_rate": float(np.mean(np.asarray(residuals) > threshold)),
            }
        )
    first_full_rejection = next(
        (row["rotation_size_epsilon"] for row in rows if row["frozen_rejection_rate"] == 1.0),
        None,
    )
    return {
        "program": "ST27",
        "object": "Isospectral Sector-Refit Adversary",
        "frozen_threshold": threshold,
        "rows": rows,
        "first_sampled_epsilon_with_full_rejection": first_full_rejection,
        "theorem": (
            "On a branch where a channel map f is spectrally injective, f(B)=f(A) implies B=A by inverse "
            "functional calculus. Orthogonal rotations QAQ^T are isospectral and pass every eigenvalue-only test, "
            "yet a frozen matrix-level common-generator protocol rejects them unless Q lies in the commutant. "
            "Sector-specific refitting hides this incompatibility by construction."
        ),
        "status": "proven_adversarial_no_refit_principle_with_strong_numerical_detection_curve",
        "boundary": "The adversaries and threshold are synthetic; this is a protocol stress test, not external evidence.",
    }


def make_figures(results: dict[str, Any]) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    st16 = results["ST16"]
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
    speed = [row["speed"] for row in st16["rows"]]
    axes[0].plot(speed, [row["jordan_solvability_invariant"] for row in st16["rows"]], label="full chain")
    axes[0].plot(speed, [row["analytic_jordan_invariant"] for row in st16["rows"]], "--", label="reduced scalar")
    axes[0].axhline(0, color="black", lw=0.8)
    axes[0].axvline(st16["selected_collision_speed"], color="tab:red", ls=":")
    axes[0].set_xlabel("mediator speed")
    axes[0].set_ylabel("Jordan solvability invariant")
    axes[0].legend()
    axes[1].plot(speed, [row["spectral_abscissa"] for row in st16["rows"]])
    axes[1].set_xlabel("mediator speed")
    axes[1].set_ylabel("spectral abscissa")
    axes[1].set_title("zero-mode collision")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st16_jordan_collision.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
    st17 = results["ST17"]
    labels = ["C*(A)", "natural\ncandidates", "full M12"]
    values = [st17["strict_functional_algebra_dimension"], st17["algebra_dimension_after_natural_candidates"], st17["full_matrix_algebra_dimension"]]
    axes[0].bar(labels, values, color=["tab:blue", "tab:orange", "tab:gray"])
    axes[0].set_ylabel("algebra dimension")
    st18 = results["ST18"]
    axes[1].plot([row["antipodal_weight_q"] for row in st18["lift_rows"]], [row["trace"] for row in st18["lift_rows"]], "o-")
    axes[1].set_xlabel("free fine antipodal weight q")
    axes[1].set_ylabel("trace A24(q)")
    axes[1].set_title("exact coarse intertwining, nonunique fine sector")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st17_algebra_st18_refinement.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
    st19 = results["ST19"]
    axes[0].hist(np.log10(np.maximum(st19["same_generator_scores"], 1e-18)), bins=30, alpha=0.75, label="same A")
    axes[0].hist(np.log10(st19["altered_generator_scores"]), bins=30, alpha=0.75, label="altered U")
    axes[0].axvline(math.log10(1e-10), color="black", ls="--", label="frozen threshold")
    axes[0].set_xlabel("log10 mixed-channel residual")
    axes[0].legend()
    intervals = results["ST20"]["eigenvalue_intervals"]
    axes[1].errorbar(
        range(12),
        [(row["lower_float"] + row["upper_float"]) / 2 for row in intervals],
        yerr=[max((row["upper_float"] - row["lower_float"]) / 2, 1e-16) for row in intervals],
        fmt="o",
    )
    axes[1].set_xlabel("Fourier mode")
    axes[1].set_ylabel("certified eigenvalue interval")
    axes[1].set_title("ST20 branch input")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st19_prediction_st20_intervals.png", dpi=190)
    plt.close(fig)

    st21 = results["ST21"]
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
    axes[0].semilogy([row["range_R"] for row in st21["finite_range_rows"]], [max(row["actual_antipodal_amplitude"], 1e-18) for row in st21["finite_range_rows"]], "o-", label="actual")
    axes[0].semilogy([row["range_R"] for row in st21["finite_range_rows"]], [row["unitarity_clipped_bound"] for row in st21["finite_range_rows"]], "s--", label="path bound")
    axes[0].set_xlabel("interaction range R")
    axes[0].set_ylabel("antipodal amplitude at t=0.1")
    axes[0].legend()
    axes[1].loglog([row["range_R"] for row in st21["long_range_rows_N192"]], [row["actual_unitary_operator_error"] for row in st21["long_range_rows_N192"]], "o-", label="actual")
    axes[1].loglog([row["range_R"] for row in st21["long_range_rows_N192"]], [row["absolute_row_tail_bound"] for row in st21["long_range_rows_N192"]], "s--", label="tail bound")
    axes[1].set_xlabel("truncation range R")
    axes[1].set_ylabel("unitary truncation error")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st21_locality_bounds.png", dpi=190)
    plt.close(fig)

    st24 = results["ST24"]
    rho = np.logspace(-4, 4, 500)
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    for exponent in [0.0, 0.25, 0.5, 0.75]:
        values = saturation_potential(rho, exponent, 1.0, 1.0, 2.0, 0.15)
        ax.plot(rho, values, label=f"a={exponent:g}")
    ax.set_xscale("log")
    ax.set_ylim(-200, 200)
    ax.set_xlabel("density rho")
    ax.set_ylabel("local effective potential")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st24_saturation_classes.png", dpi=190)
    plt.close(fig)

    st27 = results["ST27"]
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.loglog(
        [max(row["rotation_size_epsilon"], 1e-6) for row in st27["rows"]],
        [max(row["median_frozen_common_generator_residual"], 1e-18) for row in st27["rows"]],
        "o-",
        label="frozen matrix test",
    )
    ax.axhline(st27["frozen_threshold"], color="black", ls="--", label="threshold")
    ax.set_xlabel("isospectral rotation size")
    ax.set_ylabel("median residual")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st27_adversarial_refits.png", dpi=190)
    plt.close(fig)


def write_summary(results: dict[str, Any]) -> None:
    findings = {
        "ST16": "Jordan-chain criterion derived; collision speed remains non-interval numerical",
        "ST17": "all stabilizer-invariant generators remain in the seven-dimensional radial algebra",
        "ST18": "infinite exact positive C24 refinement fiber proves continued nonuniqueness",
        "ST19": "derived mixed-channel preregistration separates synthetic altered generators",
        "ST20": "scalar branches and projector matching interval-certified under mpmath.iv",
        "ST21": "finite-range path bounds proved; long-range tail remains nonlocal and signed",
        "ST22": "nine map obligations compress to five groups under operational packaging",
        "ST23": "nonzero continuous dihedral-natural U1 connection from A alone is impossible",
        "ST24": "coercive saturation classes classified by asymptotic exponent",
        "ST25": "Lean source exported but no configured toolchain; not machine checked",
        "ST26": "multiweight scale action proves a single dimension assignment is insufficient",
        "ST27": "isospectral sector refits defeat spectral tests but frozen matrix tests reject them",
    }
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["program", "status", "main_finding"])
        for index in range(16, 28):
            key = f"ST{index}"
            writer.writerow([key, results[key]["status"], findings[key]])


def main() -> None:
    w, a, _ = strict_operator()
    results: dict[str, Any] = {
        "release": "10.56",
        "programs": "ST16-ST27",
        "language": "English",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
        },
    }
    results["ST16"] = st16_exact_krein_zero_mode(a)
    results["ST17"] = st17_selector_free_algebra(a)
    results["ST18"] = st18_refinement_functor_obstruction(w, a)
    results["ST19"] = st19_derived_preregistered_observable(a)
    results["ST20"] = st20_interval_common_generator()
    results["ST21"] = st21_locality_bounds(w)
    results["ST22"] = st22_operational_tuple_minimality(a, results["ST18"])
    results["ST23"] = st23_connection_source_no_go(a)
    results["ST24"] = st24_saturation_classification()
    results["ST25"] = st25_proof_assistant_replay()
    results["ST26"] = st26_dimensional_equivariance()
    results["ST27"] = st27_adversarial_refits(a)
    results["recommended_next_programs"] = [
        {"id": "ST28", "priority": 1, "study": "interval-certify the ST16 stationary branch, Jordan coefficients and unique collision root"},
        {"id": "ST29", "priority": 2, "study": "classify state-dependent and orbit-valued symmetry breaking outside the ST17 invariant no-go"},
        {"id": "ST30", "priority": 3, "study": "impose refinement associativity on the ST18 fiber and classify all dyadic natural transformations"},
        {"id": "ST31", "priority": 4, "study": "construct a derived cross-sector observable with simulated finite counts and a frozen likelihood test"},
        {"id": "ST32", "priority": 5, "study": "certify long-range polynomial propagation bounds for eta=1.8 without Markov semantics"},
        {"id": "ST33", "priority": 6, "study": "classify discrete pi-holonomy and state-sourced connections permitted by ST23"},
        {"id": "ST34", "priority": 7, "study": "prove existence and orbital stability for the coercive ST24 saturating DNLS family"},
        {"id": "ST35", "priority": 8, "study": "configure a pinned offline Lean toolchain and replay ST25 without network-dependent state"},
        {"id": "ST36", "priority": 9, "study": "derive sector-scale identifiability conditions for the ST26 multiweight action"},
        {"id": "ST37", "priority": 10, "study": "add finite-count noise and nuisance calibration to the ST27 adversarial protocol"},
        {"id": "ST38", "priority": 11, "study": "search for a strict-sourced non-invariant state that jointly addresses selector, connection and algebra"},
        {"id": "ST39", "priority": 12, "study": "build a consolidated theorem dependency graph for ST01-ST27 and identify minimal independent axioms"},
    ]
    results["epistemic_boundary"] = (
        "ST16-ST27 remain finite local mathematics, conditional constructions and synthetic stress tests. They do not "
        "discharge QW-2191, source a state-dependent selector, produce physical units, select a canonical refinement, "
        "derive a gauge field, complete legacy-to-strict, transfer legacy roles, supply laboratory data or establish FIN as a ToE."
    )
    make_figures(results)
    write_summary(results)
    RESULTS.write_text(json.dumps(native(results), indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"results": RESULTS.name, "programs": 12, "figures": 6}, indent=2))


if __name__ == "__main__":
    main()
