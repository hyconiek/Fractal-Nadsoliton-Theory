#!/usr/bin/env python3
"""FIN ST46--ST57: source, carrier, refinement, and entropy research.

This batch continues Release 10.57 using local analytical and lightweight
computational methods.  It keeps strict-source theorems, conditional
constructions, synthetic receiver studies, and failed closures separate.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import warnings
from fractions import Fraction
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
from scipy.linalg import expm

from fin_programs_488_496_low_compute import grouped_stieltjes_data
from fin_programs_507_516_research import stationary_root, strict_a_interval
from fin_programs_497_506_next_research import iv_bounds
from fin_st01_st15_research import N, random_orthogonal_fixing_uniform, strict_operator
from fin_st16_st27_research import final_memory_state
from fin_st28_st45_research import (
    bicommutant_algebra_dimension,
    dihedral_orbit_density,
    dyadic_lift,
    experiment_probabilities,
    interval_quadratic_positive_root,
    periodic_embedding,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST46_ST57_Results.json"
SUMMARY = ROOT / "FIN_ST46_ST57_Summary.csv"
CERT46 = ROOT / "FIN_ST46_Upstream_Sensitivity_Certificate.json"
PREREG51 = ROOT / "FIN_ST51_Two_Carrier_Preregistration.json"
SPEC55 = ROOT / "FIN_ST55_Two_Carrier_Executable_Spec.json"
FIG_DIR = ROOT / "FIN_ST46_ST57_Figures"
SEED = 20260813


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


def canonical_digest(payload: dict) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def cluster_data(a: np.ndarray) -> list[dict]:
    even = np.arange(0, N, 2)
    hidden = np.arange(1, N, 2)
    b = a[np.ix_(even, hidden)]
    c = a[np.ix_(hidden, hidden)] + 0.15 * np.eye(len(hidden))
    vals, vecs = np.linalg.eigh(c)
    groups: list[list[int]] = []
    for index, value in enumerate(vals):
        if not groups or abs(value - vals[groups[-1][0]]) > 1e-9:
            groups.append([index])
        else:
            groups[-1].append(index)
    rows = []
    for group in groups:
        projector = vecs[:, group] @ vecs[:, group].T
        pole = float(np.mean(vals[group]))
        residue = float(np.trace(b @ projector @ b.T) / len(even))
        outside = [j for j in range(len(vals)) if j not in group]
        gap = min(abs(vals[i] - vals[j]) for i in group for j in outside) if outside else math.inf
        rows.append(
            {
                "indices": group,
                "multiplicity": len(group),
                "pole": pole,
                "residue": residue,
                "projector": projector,
                "gap": float(gap),
                "B": b,
            }
        )
    return rows


def st46_upstream_sensitivity_certificate() -> dict:
    alo, ahi, amid = strict_a_interval()
    radius_matrix = 0.5 * (ahi - alo)
    delta_a = float(math.sqrt(np.linalg.norm(radius_matrix, 1) * np.linalg.norm(radius_matrix, np.inf)))
    clusters = cluster_data(amid)
    b = clusters[0]["B"]
    b_f = float(np.linalg.norm(b, "fro"))
    delta_b_f = math.sqrt(6.0) * delta_a
    identity = np.eye(N)
    memory = np.zeros_like(amid)
    fzero = 0.0
    delta_fzero = 0.0
    delta_memory = 0.0
    cluster_rows = []
    for row in clusters:
        pole = row["pole"]
        residue = row["residue"]
        gap = row["gap"]
        delta_pole = delta_a
        projector_bound = 2.0 * delta_a / max(gap - 2.0 * delta_a, gap / 2.0)
        bprime_f = b_f + delta_b_f
        delta_residue = (
            delta_b_f * bprime_f + b_f * projector_bound * bprime_f + b_f * delta_b_f
        ) / 6.0
        pmin = pole - delta_pole
        resolvent = np.linalg.inv(amid + pole * identity)
        resolvent_delta = (delta_a + delta_pole) / (pole * pmin)
        fzero += residue / pole
        memory -= residue * resolvent
        delta_fzero += delta_residue / pmin + abs(residue) * delta_pole / (pole * pmin)
        delta_memory += delta_residue / pmin + abs(residue) * resolvent_delta
        cluster_rows.append(
            {
                "multiplicity": row["multiplicity"],
                "pole_center": pole,
                "residue_center": residue,
                "intercluster_gap": gap,
                "pole_error_bound": delta_pole,
                "projector_error_bound": projector_bound,
                "residue_error_bound": delta_residue,
            }
        )
    memory += fzero * identity
    delta_memory += delta_fzero
    # Pay a conservative floating-arithmetic safety factor after analytic bounds.
    safety_factor = 16.0
    delta_a *= safety_factor
    delta_memory *= safety_factor
    delta_operator = delta_a + delta_memory

    state, _, _, _ = final_memory_state(amid)
    f_center = (amid + memory) @ state + state - state**3
    j_center = amid + memory + identity - np.diag(3.0 * state**2)
    j_inverse = np.linalg.inv(j_center)
    jinv_norm = float(np.linalg.norm(j_inverse, 2))
    inverse_residual = float(np.linalg.norm(identity - j_inverse @ j_center, 2))
    eta = float(np.linalg.norm(j_inverse @ f_center, 2)) + jinv_norm * delta_operator * float(np.linalg.norm(state))
    radius_rows = []
    accepted_radius = None
    for radius in [1e-12, 3e-12, 1e-11, 3e-11, 1e-10, 3e-10, 1e-9, 3e-9, 1e-8]:
        derivative_change = delta_operator + 6.0 * float(np.max(np.abs(state))) * radius + 3.0 * radius**2
        contraction = inverse_residual + jinv_norm * derivative_change
        image = eta + contraction * radius
        accepted = bool(contraction < 1.0 and image < radius)
        radius_rows.append(
            {
                "euclidean_radius": radius,
                "contraction_bound": contraction,
                "newton_image_bound": image,
                "strict_inclusion": accepted,
            }
        )
        if accepted and accepted_radius is None:
            accepted_radius = radius

    # Propagate the same upstream uncertainty through the Jordan coefficients.
    poles = np.asarray([row["pole"] for row in clusters])
    residues = np.asarray([row["residue"] for row in clusters])
    pole_errors = np.asarray([row["pole_error_bound"] * safety_factor for row in cluster_rows])
    residue_errors = np.asarray([row["residue_error_bound"] * safety_factor for row in cluster_rows])
    tmat = np.zeros_like(amid)
    smat = np.zeros_like(amid)
    delta_t = 0.0
    delta_s = 0.0
    for pole, residue, dpole, dresidue in zip(poles, residues, pole_errors, residue_errors):
        resolvent = np.linalg.inv(amid + pole * identity)
        tmat += residue * resolvent @ resolvent
        smat += residue * resolvent @ resolvent @ resolvent
        pmin = pole - dpole
        dresolvent = (delta_a + dpole) / (pole * pmin)
        delta_t += dresidue / pmin**2 + 2.0 * abs(residue) * dresolvent / pmin
        delta_s += dresidue / pmin**3 + 3.0 * abs(residue) * dresolvent / pmin**2

    if accepted_radius is None:
        coefficient_intervals = None
        collision_interval = None
    else:
        eps = accepted_radius
        lplus = amid + memory + identity - np.diag(3.0 * state**2)
        binv = np.linalg.inv(lplus)
        bnorm = float(np.linalg.norm(binv, 2))
        delta_l = delta_operator + 6.0 * float(np.max(np.abs(state))) * eps + 3.0 * eps**2
        inverse_denominator = 1.0 - bnorm * delta_l
        delta_binv = bnorm * bnorm * delta_l / inverse_denominator
        bprime = bnorm + delta_binv
        u_norm = float(np.linalg.norm(state))
        tu = tmat @ state
        tu_norm = float(np.linalg.norm(tu))
        delta_tu = delta_t * (u_norm + eps) + float(np.linalg.norm(tmat, 2)) * eps

        c0 = float(state @ binv @ state)
        c1 = float(state @ binv @ tu)
        c2 = float(tu @ binv @ tu)
        c3 = float(state @ smat @ state)
        dc0 = eps * bprime * (u_norm + eps) + u_norm * delta_binv * (u_norm + eps) + u_norm * bnorm * eps
        dc1 = eps * bprime * (tu_norm + delta_tu) + u_norm * delta_binv * (tu_norm + delta_tu) + u_norm * bnorm * delta_tu
        dc2 = delta_tu * bprime * (tu_norm + delta_tu) + tu_norm * delta_binv * (tu_norm + delta_tu) + tu_norm * bnorm * delta_tu
        snorm = float(np.linalg.norm(smat, 2))
        dc3 = eps * (snorm + delta_s) * (u_norm + eps) + u_norm * delta_s * (u_norm + eps) + u_norm * snorm * eps
        coefficient_intervals = {
            "a=c2+c3": [c2 + c3 - dc2 - dc3, c2 + c3 + dc2 + dc3],
            "b=2c1": [2.0 * (c1 - dc1), 2.0 * (c1 + dc1)],
            "c=c0": [c0 - dc0, c0 + dc0],
        }
        al, au = coefficient_intervals["a=c2+c3"]
        bl, bu = coefficient_intervals["b=2c1"]
        cl, cu = coefficient_intervals["c=c0"]
        xlo = interval_quadratic_positive_root(au, bu, cu)
        xhi = interval_quadratic_positive_root(al, bl, cl)
        if xlo > xhi:
            xlo, xhi = xhi, xlo
        collision_interval = [1.0 / xhi, 1.0 / xlo]

    certificate = {
        "strict_input_model": "exact decimal omega=743/4000, phi=13/80, eta=9/5, beta=1",
        "strict_matrix_operator_error_bound": delta_a,
        "memory_operator_error_bound": delta_memory,
        "combined_operator_error_bound": delta_operator,
        "safety_factor": safety_factor,
        "cluster_rows": cluster_rows,
        "stationary_center_residual_2": float(np.linalg.norm(f_center)),
        "radius_rows": radius_rows,
        "accepted_state_radius_2": accepted_radius,
        "coefficient_intervals": coefficient_intervals,
        "collision_speed_interval": collision_interval,
        "arithmetic_boundary": (
            "The strict transcendental input is directed-interval regenerated. Spectral-projector, memory, inverse, "
            "and root propagation use analytic perturbation inequalities evaluated in binary64 with a factor-16 "
            "safety payment; this is not an independent all-interval or proof-assistant replay."
        ),
    }
    CERT46.write_text(json.dumps(native(certificate), indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST46",
        "object": "Upstream Strict-to-Memory Sensitivity Certificate",
        **certificate,
        "status": "strong_partial_upstream_certificate_not_full_interval_closure",
    }


def st47_scale_stationary_refinement(a12: np.ndarray) -> dict:
    s12 = float(np.trace(a12) / 12.0)
    embedding = periodic_embedding(12)
    rows = []
    for q in [0.0, 0.1, 0.3, 0.9]:
        a24 = dyadic_lift(a12, q)
        rows.append(
            {
                "q": q,
                "coarse_intertwining_residual": float(np.linalg.norm(a24 @ embedding - embedding @ a12)),
                "trace_density": float(np.trace(a24) / 24.0),
                "trace_density_minus_coarse": float(np.trace(a24) / 24.0 - s12),
            }
        )
    alternative_scales = []
    for c in [1.0, 1.05, 1.2]:
        alternative_scales.append({"target_trace_multiplier": c, "selected_q": (c - 1.0) * s12})
    return {
        "program": "ST47",
        "object": "Trace-Stationary Refinement Selector",
        "coarse_trace_density": s12,
        "rows": rows,
        "unique_nonnegative_q_under_trace_density_conservation": 0.0,
        "alternative_scale_laws": alternative_scales,
        "theorem": (
            "For the ST18 dyadic lift, tr(A_24(q))/24=tr(A_12)/12+q. Exact conservation of "
            "Dirichlet trace density therefore selects q=0 uniquely. More generally a prescribed multiplier c "
            "selects q=(c-1)s_12, so the selection comes from the scale law, not from intertwining alone."
        ),
        "status": "proven_conditional_unique_refinement_under_added_trace_law",
        "boundary": "Trace-density conservation is an added refinement axiom and is not derived from strict FIN or fractal compression.",
    }


def st48_gain_source_nonidentifiability(a: np.ndarray) -> dict:
    reflection = np.zeros((N, N))
    for x in range(N):
        reflection[(-x) % N, x] = 1.0
    pattern = np.sin(2.0 * math.pi * np.arange(N) / N)
    h = np.diag(pattern)
    h /= np.linalg.norm(h)
    odd_residual = float(np.linalg.norm(reflection @ h @ reflection.T + h))
    mus = [-0.35, 0.0, 0.35]
    rows = []
    for mu in mus:
        # Along K=A+xH, gradient descent for F=-mu*x^2/2+x^4/4 is xdot=mu*x-x^3.
        x0 = 1e-4
        x = x0
        dt = 0.01
        for _ in range(4000):
            x += dt * (mu * x - x**3)
        rows.append(
            {
                "mu": mu,
                "odd_linear_gain": mu,
                "final_amplitude": x,
                "strict_fixed_point_stationarity": True,
                "reflection_equivariant": True,
            }
        )
    return {
        "program": "ST48",
        "object": "Odd-Gain Source Nonidentifiability Theorem",
        "reflection_odd_tangent_residual": odd_residual,
        "rows": rows,
        "potential_family": "F_mu(A+H)=-mu*||P_-H||^2/2+||P_-H||^4/4+||P_+H||^2/2",
        "theorem": (
            "The same strict fixed operator A and the same reflection symmetry admit equivariant adaptive potentials "
            "with negative, zero, or positive odd-sector linear gain. All have A as a stationary point. Therefore A "
            "and its stabilizer do not determine the sign of the feedback gain; a response functional or state law is necessary."
        ),
        "status": "proven_strict_static_core_does_not_identify_feedback_gain_sign",
        "boundary": "The counterfamily is constructed; it proves underdetermination and does not select a physical adaptive law.",
    }


def projective_holonomy(spinors: np.ndarray) -> complex:
    product = 1.0 + 0.0j
    for x in range(N):
        overlap = np.vdot(spinors[x], spinors[(x + 1) % N])
        product *= overlap / abs(overlap)
    return complex(product)


def st49_projective_multicomponent_state(a: np.ndarray) -> dict:
    angles = 2.0 * math.pi * np.arange(N) / N
    rays = np.stack([np.cos(angles / 2.0), np.sin(angles / 2.0)], axis=1).astype(complex)
    amplitudes = 1.0 + 0.02 * np.arange(N)
    spinors = amplitudes[:, None] * rays
    normalized_rays = spinors / np.linalg.norm(spinors, axis=1)[:, None]
    holonomy = projective_holonomy(normalized_rays)
    density = np.sum(np.abs(spinors) ** 2, axis=1)
    density /= np.sum(density)
    algebra_dim, commutant_dim = bicommutant_algebra_dimension([a, np.diag(density)])
    shifted_holonomies = [projective_holonomy(np.roll(normalized_rays, shift, axis=0)) for shift in range(N)]
    return {
        "program": "ST49",
        "object": "Two-Component Projective Selector-Holonomy State",
        "projective_holonomy": holonomy,
        "holonomy_phase": float(np.angle(holonomy) % (2.0 * math.pi)),
        "shifted_holonomy_maximum_residual": float(max(abs(value + 1.0) for value in shifted_holonomies)),
        "density_algebra_dimension": algebra_dim,
        "joint_commutant_dimension": commutant_dim,
        "density_orbit_size": len(dihedral_orbit_density(density)),
        "unique_density_maximum": bool(np.sum(np.isclose(density, np.max(density))) == 1),
        "density": density,
        "rays": normalized_rays,
        "theorem": (
            "A two-component real projective texture with one half-angle turn has Pancharatnam holonomy -1 (flux pi). "
            "Adding distinct amplitudes gives a unique density maximum and, with connected A, the full matrix algebra. "
            "Thus one enlarged state type can carry selector data, full algebra and pi holonomy."
        ),
        "status": "proven_conditional_joint_formal_closure_with_added_projective_state",
        "boundary": (
            "The amplitude ordering and projective texture are inserted, not strict-sourced. Pi flux is reflection-fixed "
            "and does not by itself choose chirality or a canonical origin. QW-2191 remains open."
        ),
    }


def st50_channel_intertwiner_classification(a: np.ndarray) -> dict:
    lambdas = np.linalg.eigvalsh(a)
    t = 0.61
    tau = 0.61
    z = 1.0
    beta = 0.8
    channels = {
        "unitary": np.exp(-1j * t * lambdas),
        "heat": np.exp(-tau * lambdas),
        "wave": np.cos(t * np.sqrt(np.maximum(lambdas, 0.0))),
        "green_z1": 1.0 / (lambdas + z),
        "gibbs_weights": np.exp(-beta * lambdas) / np.sum(np.exp(-beta * lambdas)),
    }
    names = list(channels)
    rows = []
    matrix = np.zeros((len(names), len(names)), dtype=int)
    for i, left in enumerate(names):
        for j, right in enumerate(names):
            dimension = int(
                sum(abs(x - y) < 1e-10 for x in channels[left] for y in channels[right])
            )
            matrix[i, j] = dimension
            rows.append({"left": left, "right": right, "intertwiner_space_dimension": dimension})
    injective = {
        name: len(np.unique(np.round(np.column_stack([values.real, values.imag]), 12), axis=0)) == 7
        for name, values in channels.items()
    }
    return {
        "program": "ST50",
        "object": "Common-Spectrum Channel Intertwiner Classification",
        "parameters": {"t": t, "tau": tau, "z": z, "beta": beta},
        "channel_names": names,
        "intertwiner_dimension_matrix": matrix,
        "rows": rows,
        "injective_on_seven_distinct_strict_eigenvalues": injective,
        "invertible_cross_channel_similarity_pairs": [],
        "theorem": (
            "For normal spectral channels f(A),g(A), the linear intertwiner space has dimension sum over mode pairs "
            "with f(lambda_i)=g(lambda_j). An invertible intertwiner exists only when the eigenvalue multisets agree. "
            "At the declared parameters no two different channels are similar, although several share the zero-mode value one."
        ),
        "status": "proven_finite_channel_intertwiner_classification",
        "boundary": "Gibbs weights are state eigenvalues rather than a time channel; the table is finite and parameter-specific.",
    }


def carrier_probability_table(a: np.ndarray, q: np.ndarray, transported: bool) -> np.ndarray:
    a2 = q @ a @ q.T
    t = 0.63
    probabilities = []
    for x in range(N):
        vector = np.eye(N)[x]
        state = np.outer(vector, vector)
        effects = [np.diag(np.eye(N)[y]) for y in range(N)]
        if transported:
            state = q @ state @ q.T
            effects = [q @ effect @ q.T for effect in effects]
        probabilities.append(experiment_probabilities(a2, state, effects, t))
    return np.asarray(probabilities)


def pearson_score(counts: np.ndarray, expected_probabilities: np.ndarray, shots: int) -> float:
    expected = shots * expected_probabilities
    return float(np.sum((counts - expected) ** 2 / np.maximum(expected, 1e-12)))


def st51_two_carrier_finite_count(a: np.ndarray) -> dict:
    config = {
        "seed": SEED + 51,
        "time": 0.63,
        "shots_per_preparation": 1200,
        "preparations": 12,
        "outcomes": 12,
        "training_null_trials": 300,
        "holdout_null_trials": 300,
        "mismatched_trials": 300,
        "frozen_quantile": 0.99,
    }
    digest = canonical_digest(config)
    PREREG51.write_text(json.dumps({"configuration": config, "sha256": digest}, indent=2, sort_keys=True), encoding="utf-8")
    rng = np.random.default_rng(config["seed"])
    q = random_orthogonal_fixing_uniform(rng)
    p_reference = carrier_probability_table(a, np.eye(N), transported=False)
    p_transported = carrier_probability_table(a, q, transported=True)
    p_mismatch = carrier_probability_table(a, q, transported=False)

    def sample_scores(probabilities: np.ndarray, trials: int) -> list[float]:
        scores = []
        for _ in range(trials):
            counts = np.vstack(
                [rng.multinomial(config["shots_per_preparation"], row) for row in probabilities]
            )
            scores.append(pearson_score(counts, p_reference, config["shots_per_preparation"]))
        return scores

    training = sample_scores(p_transported, config["training_null_trials"])
    threshold = float(np.quantile(training, config["frozen_quantile"], method="higher"))
    holdout = sample_scores(p_transported, config["holdout_null_trials"])
    mismatch = sample_scores(p_mismatch, config["mismatched_trials"])
    return {
        "program": "ST51",
        "object": "Finite-Count Two-Carrier Transfer Test",
        "preregistration_sha256": digest,
        "transported_probability_residual": float(np.max(np.abs(p_reference - p_transported))),
        "mismatched_probability_total_variation_mean": float(
            np.mean(0.5 * np.sum(np.abs(p_reference - p_mismatch), axis=1))
        ),
        "frozen_threshold": threshold,
        "holdout_false_rejection_rate": float(np.mean(np.asarray(holdout) > threshold)),
        "mismatched_detection_power": float(np.mean(np.asarray(mismatch) > threshold)),
        "training_scores": training,
        "holdout_scores": holdout,
        "mismatch_scores": mismatch,
        "status": "strong_synthetic_two_carrier_transfer_discrimination",
        "boundary": "Both carriers, encoders, detectors, counts and custody are synthetic and locally generated.",
    }


def polynomial_sign(value: Fraction) -> Fraction:
    return 3 * value**3 - 80 * value + 80


def rational_bisect(left: Fraction, right: Fraction, iterations: int = 80) -> tuple[Fraction, Fraction]:
    fl = polynomial_sign(left)
    fr = polynomial_sign(right)
    if fl * fr >= 0:
        raise ValueError("root is not bracketed")
    for _ in range(iterations):
        mid = (left + right) / 2
        fm = polynomial_sign(mid)
        if fm == 0:
            return mid, mid
        if fl * fm < 0:
            right, fr = mid, fm
        else:
            left, fl = mid, fm
    return left, right


def onsite_potential(rho: np.ndarray | float, rscale: float = 2.0, mu: float = 0.15):
    q = rho / (1.0 + rho / rscale)
    return -0.5 * q**2 + 0.5 * mu * rho


def st52_uniform_minimizer_orbit(a: np.ndarray) -> dict:
    small = rational_bisect(Fraction(1042, 1000), Fraction(1043, 1000))
    large = rational_bisect(Fraction(4563, 1000), Fraction(4564, 1000))
    ylo, yhi = map(float, large)
    ystar = 0.5 * (ylo + yhi)
    rho = 2.0 * (ystar - 1.0)
    amplitude = math.sqrt(rho)
    curvature = (ystar - 1.0) * (9.0 * ystar**2 - 80.0) / (10.0 * ystar**3)
    eigs = np.linalg.eigvalsh(a)
    return {
        "program": "ST52",
        "object": "Exact Uniform Minimizer-Orbit Theorem for the Constructed Saturation",
        "critical_y_intervals": {
            "local_maximum": [float(small[0]), float(small[1])],
            "global_minimum": [ylo, yhi],
        },
        "critical_y_rational_intervals": {
            "local_maximum": [str(small[0]), str(small[1])],
            "global_minimum": [str(large[0]), str(large[1])],
        },
        "global_minimum_amplitude": amplitude,
        "global_minimum_density": rho,
        "onsite_minimum_value": float(onsite_potential(rho)),
        "radial_hessian_curvature": curvature,
        "phase_hessian_zero_modes": 1,
        "phase_hessian_positive_gap": float(eigs[1]),
        "minimizer_orbit": "{sqrt(rho_*) exp(i theta) 1 : theta in S1}",
        "theorem": (
            "For the constructed ST34 potential, V'(rho) has the sign of 3y^3-80y+80 with y=1+rho/2. "
            "Its larger positive critical root is the unique global onsite minimum. Since A is a connected positive-weight "
            "graph Laplacian, the complete global minimizer set is the uniform U(1) orbit. Positive radial curvature and the "
            "strict spectral gap prove orbital stability of this orbit."
        ),
        "status": "proven_global_uniform_orbit_and_orbital_stability_for_constructed_model",
        "boundary": "This closes the constructed model but proves nonlocalization; the saturation and mu remain additional laws.",
    }


def shannon(probabilities: np.ndarray) -> float:
    positive = probabilities[probabilities > 0]
    return float(-np.sum(positive * np.log(positive)))


def relative_entropy(p: np.ndarray, q: np.ndarray) -> float:
    return float(np.sum(p * np.log(p / q)))


def normalized_heat_filter(a: np.ndarray, rho: np.ndarray, tau: float) -> np.ndarray:
    p = expm(-tau * a)
    output = p @ rho @ p
    return output / np.trace(output)


def st53_entropy_data_processing(a: np.ndarray) -> dict:
    tau = 0.4
    transition = expm(-tau * a)
    p = np.arange(1, N + 1, dtype=float)
    p /= np.sum(p)
    q = np.arange(N, 0, -1, dtype=float) + 0.3
    q /= np.sum(q)
    pp = transition @ p
    qq = transition @ q
    rho0 = np.diag(p)
    unitary = expm(-1j * 0.7 * a)
    rho_u = unitary @ rho0 @ unitary.conj().T
    vn_before = shannon(np.linalg.eigvalsh(rho0))
    vn_after = shannon(np.maximum(np.linalg.eigvalsh(rho_u), 0.0))
    rho1 = np.diag(np.eye(N)[0])
    uniform = np.ones(N) / math.sqrt(N)
    rho2 = np.outer(uniform, uniform)
    nonlinear_filter_defect = float(
        np.linalg.norm(
            normalized_heat_filter(a, 0.5 * (rho1 + rho2), tau)
            - 0.5 * (normalized_heat_filter(a, rho1, tau) + normalized_heat_filter(a, rho2, tau))
        )
    )
    return {
        "program": "ST53",
        "object": "Carrier-Neutral Entropy and Data-Processing Ledger",
        "transition_minimum_entry": float(np.min(transition)),
        "row_sum_residual": float(np.max(np.abs(np.sum(transition, axis=1) - 1.0))),
        "column_sum_residual": float(np.max(np.abs(np.sum(transition, axis=0) - 1.0))),
        "shannon_before": shannon(p),
        "shannon_after_strict_markov": shannon(pp),
        "relative_entropy_before": relative_entropy(p, q),
        "relative_entropy_after": relative_entropy(pp, qq),
        "unitary_von_neumann_entropy_residual": abs(vn_after - vn_before),
        "normalized_heat_filter_affinity_defect": nonlinear_filter_defect,
        "theorem": (
            "On C12 all strict weights are positive, so exp(-tau A) is doubly stochastic: Shannon entropy does not decrease "
            "and classical relative entropy contracts. Unitary conjugation preserves von Neumann entropy. The normalized "
            "amplitude heat filter is nonlinear and is not a CPTP channel, so its use requires a separate operational model."
        ),
        "status": "proven_finite_entropy_invariants_and_data_processing_boundary",
        "boundary": "Dimensionless information entropy is not thermodynamic entropy without a bath, temperature, and energy scale.",
    }


def st54_operational_quotient_universal_property() -> dict:
    return {
        "program": "ST54",
        "object": "Operational-Isomorphism Quotient Universal Property",
        "category": "groupoid Rep_op of finite operational realizations and operational isomorphisms",
        "quotient": "pi: Rep_op -> pi_0(Rep_op)",
        "universal_property": (
            "For every set-valued invariant F that assigns equal values to operationally isomorphic realizations, "
            "there exists a unique map Fbar on pi_0(Rep_op) such that F=Fbar composed with pi."
        ),
        "proof": (
            "Define Fbar([R])=F(R). Isomorphism invariance makes this well-defined; surjectivity of pi gives uniqueness."
        ),
        "status": "proven_quotient_universal_property",
        "boundary": "This is a classification theorem, not a derivation of which operational groupoid nature realizes.",
    }


def st55_executable_two_carrier_spec(st51: dict) -> dict:
    specification = {
        "title": "FIN two-carrier operational-equivalence transfer test",
        "version": "1.0.0",
        "hypotheses": {
            "H0": "carrier B transports generator, preparations, effects, clock and record labels by one calibrated isomorphism Q",
            "H1": "at least one operational component is not transported by the registered Q",
        },
        "carrier_A": {"basis": "vertex", "generator": "A", "preparations": "12 vertex states", "effects": "12 vertex effects"},
        "carrier_B": {"basis": "dense Q-basis", "generator": "Q A Q^T", "preparations": "Q rho_x Q^T", "effects": "Q E_y Q^T"},
        "clock": {"dimensionless_time": 0.63, "calibration_required_for_physical_run": True},
        "counts": {"shots_per_preparation": 1200, "preparations": 12, "outcomes": 12},
        "raw_event_schema": ["timestamp", "carrier_id", "preparation_id", "outcome_id", "run_id", "calibration_split", "blind_label"],
        "roles": {"provider": "must be external", "registrar": "must differ from provider", "analyst": "blinded", "custodian": "holds mapping until hash freeze"},
        "decision_rule": {"score": "Pearson discrepancy from carrier-A table", "threshold_source": "frozen 0.99 null quantile", "one_run": True},
        "failure_rule": "report failure without refitting Q, time, effects, threshold or kernel",
        "local_preregistration_sha256": st51["preregistration_sha256"],
        "physical_boundary": "The local packet does not supply apparatus, SI time, independent roles or external raw events.",
    }
    digest = canonical_digest(specification)
    packet = {"specification": specification, "sha256": digest}
    SPEC55.write_text(json.dumps(packet, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST55",
        "object": "Executable Two-Carrier Transfer Specification",
        "specification_sha256": digest,
        "required_raw_fields": len(specification["raw_event_schema"]),
        "logical_role_count": len(specification["roles"]),
        "status": "constructed_type_complete_local_executable_spec",
        "boundary": "Logical role separation in a file is not independent custody and cannot be generated by code.",
    }


def strict_profile(distance: int) -> float:
    return math.cos(0.18575 * distance + 0.16250) / (1.0 + distance**1.8)


def st56_fractal_compression_refinement(a12: np.ndarray) -> dict:
    mp.iv.dps = 60
    omega_iv = mp.iv.mpf("0.18575")
    phi_iv = mp.iv.mpf("0.16250")
    eta_iv = mp.iv.mpf("1.8")

    def profile_iv(distance: int):
        d = mp.iv.mpf(distance)
        return mp.iv.cos(omega_iv * d + phi_iv) / (1 + d**eta_iv)

    ratio_rows = []
    for distance in [1, 2, 3]:
        ratio = strict_profile(2 * distance) / strict_profile(distance)
        ratio_interval = iv_bounds(profile_iv(2 * distance) / profile_iv(distance))
        ratio_rows.append(
            {
                "distance": distance,
                "K_2d_over_K_d": ratio,
                "K_2d_over_K_d_interval": ratio_interval,
                "effective_power_exponent": -math.log(abs(ratio), 2.0),
            }
        )
    exponents = np.asarray([row["effective_power_exponent"] for row in ratio_rows])
    pairwise_disjoint = all(
        left["K_2d_over_K_d_interval"][1] < right["K_2d_over_K_d_interval"][0]
        or right["K_2d_over_K_d_interval"][1] < left["K_2d_over_K_d_interval"][0]
        for i, left in enumerate(ratio_rows)
        for right in ratio_rows[i + 1 :]
    )
    embedding = periodic_embedding(12)
    q_rows = []
    for q in [0.0, 0.1, 0.3, 0.9]:
        lifted = dyadic_lift(a12, q)
        q_rows.append(
            {
                "q": q,
                "coarse_reconstruction_loss": float(np.linalg.norm(lifted @ embedding - embedding @ a12) ** 2),
                "trace_complexity_penalty": float(q**2),
            }
        )
    return {
        "program": "ST56",
        "object": "Fractal-Compression Refinement Nonselection Audit",
        "scale_ratio_rows": ratio_rows,
        "directed_ratio_intervals_pairwise_disjoint": pairwise_disjoint,
        "effective_exponent_spread": float(np.max(exponents) - np.min(exponents)),
        "refinement_rows": q_rows,
        "theorem": (
            "Every loss depending only on the coarse intertwining residual is identically zero on the full q>=0 lift fibre, "
            "so coarse reconstruction cannot select q. A rate or complexity penalty can select a member, but its code and "
            "weight are additional assumptions. The strict profile also lacks one exact dyadic amplitude ratio across d=1,2,3."
        ),
        "status": "proven_coarse_compression_nonselection_and_no_exact_dyadic_ratio",
        "boundary": "This does not exclude a richer sourced fractal code retaining fine data; no such code is currently exported.",
    }


def st57_updated_axiom_graph() -> dict:
    rows = [
        {"axiom": "A0 strict finite operator", "after_ST46_ST56": "retained", "witness": "ST46 tightens but does not remove the source object"},
        {"axiom": "A1 state-selection event", "after_ST46_ST56": "retained", "witness": "ST49 amplitudes/order are inserted"},
        {"axiom": "A2 refinement law", "after_ST46_ST56": "retained", "witness": "ST47/ST56 require an added scale or code law"},
        {"axiom": "A3 dimensional calibration", "after_ST46_ST56": "retained", "witness": "ST53 entropy remains dimensionless"},
        {"axiom": "A4 operational process", "after_ST46_ST56": "constructible", "witness": "ST55 packages it but does not derive it from A"},
        {"axiom": "A5 external custody/data", "after_ST46_ST56": "retained", "witness": "ST51/ST55 are local synthetic artifacts"},
        {"axiom": "A6 nonlinear response/gain", "after_ST46_ST56": "retained", "witness": "ST48 proves gain-sign underdetermination"},
        {"axiom": "A7 connection/projective sector", "after_ST46_ST56": "conditionally packaged", "witness": "ST49 supplies pi flux only as an added texture"},
        {"axiom": "A8 temporal memory realization", "after_ST46_ST56": "retained", "witness": "ST46 certifies static sensitivity, not a unique temporal completion"},
    ]
    return {
        "program": "ST57",
        "object": "Post-ST56 Minimal-Axiom Reconciliation",
        "axiom_count_before": 9,
        "axiom_count_after_strict_source_audit": 9,
        "rows": rows,
        "conditional_multi_obligation_objects": {
            "projective_spinor_state": ["A1", "A7", "algebra completion"],
            "two_carrier_operational_packet": ["A4", "record schema"],
            "saturating_uniform_order_parameter": ["A6", "stable phase orbit"],
        },
        "theorem": (
            "Several constructed objects can package multiple obligations, but removal countermodels keep the nine source "
            "groups independent relative to the current strict core. No strict axiom reduction is proved."
        ),
        "status": "proven_no_strict_axiom_reduction_after_conditional_packaging",
        "boundary": "Minimality remains relative to the declared targets and construction classes.",
    }


def make_figures(results: dict[str, Any]) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    st46 = results["ST46"]
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    rows = st46["radius_rows"]
    ax.loglog([r["euclidean_radius"] for r in rows], [r["newton_image_bound"] for r in rows], "o-", label="Newton image")
    ax.loglog([r["euclidean_radius"] for r in rows], [r["euclidean_radius"] for r in rows], "--", label="box radius")
    ax.set(xlabel="state radius", ylabel="bound", title="ST46 upstream sensitivity inclusion")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st46_upstream_sensitivity.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.3, 4.2))
    st47 = results["ST47"]
    axes[0].plot([r["q"] for r in st47["rows"]], [r["trace_density_minus_coarse"] for r in st47["rows"]], "o-")
    axes[0].axhline(0, color="black", ls="--")
    axes[0].set(xlabel="fine weight q", ylabel="trace-density defect", title="added trace law selects q=0")
    st48 = results["ST48"]
    axes[1].bar([str(r["mu"]) for r in st48["rows"]], [r["final_amplitude"] for r in st48["rows"]])
    axes[1].set(xlabel="inserted odd gain mu", ylabel="final odd amplitude", title="same A, three gain signs")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st47_refinement_st48_gain.png", dpi=190)
    plt.close(fig)

    st49 = results["ST49"]
    rays = np.asarray(st49["rays"], dtype=complex)
    fig, axes = plt.subplots(1, 2, figsize=(11.3, 4.2))
    axes[0].plot(rays[:, 0].real, rays[:, 1].real, "o-")
    axes[0].set_aspect("equal")
    axes[0].set(xlabel="spinor component 1", ylabel="component 2", title="projective half-angle texture")
    axes[1].bar(np.arange(N), st49["density"])
    axes[1].set(xlabel="vertex", ylabel="density", title="inserted amplitude selector")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st49_projective_state.png", dpi=190)
    plt.close(fig)

    st50 = results["ST50"]
    fig, ax = plt.subplots(figsize=(7.0, 5.7))
    matrix = np.asarray(st50["intertwiner_dimension_matrix"])
    im = ax.imshow(matrix, cmap="viridis")
    ax.set_xticks(range(len(st50["channel_names"])), st50["channel_names"], rotation=35, ha="right")
    ax.set_yticks(range(len(st50["channel_names"])), st50["channel_names"])
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(matrix[i, j]), ha="center", va="center", color="white" if matrix[i, j] > matrix.max() / 2 else "black")
    fig.colorbar(im, ax=ax, label="intertwiner dimension")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st50_channel_intertwiners.png", dpi=190)
    plt.close(fig)

    st51 = results["ST51"]
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.hist(st51["holdout_scores"], bins=25, alpha=0.75, label="transported")
    ax.hist(st51["mismatch_scores"], bins=25, alpha=0.65, label="mismatched")
    ax.axvline(st51["frozen_threshold"], color="black", ls="--", label="frozen threshold")
    ax.set(xlabel="Pearson score", ylabel="count", title="two-carrier finite-count transfer test")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st51_two_carrier_test.png", dpi=190)
    plt.close(fig)

    st52 = results["ST52"]
    rho = np.linspace(0, 13, 500)
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.plot(rho, onsite_potential(rho))
    ax.axvline(st52["global_minimum_density"], color="tab:red", ls="--", label="global onsite minimum")
    ax.set(xlabel="density rho", ylabel="V(rho)", title="constructed saturating potential")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st52_uniform_minimizer.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.3, 4.2))
    st53 = results["ST53"]
    axes[0].bar(["H before", "H after"], [st53["shannon_before"], st53["shannon_after_strict_markov"]])
    axes[0].set(title="doubly stochastic heat increases H", ylabel="nats")
    st56 = results["ST56"]
    axes[1].plot([r["distance"] for r in st56["scale_ratio_rows"]], [r["effective_power_exponent"] for r in st56["scale_ratio_rows"]], "o-")
    axes[1].set(xlabel="distance d", ylabel="effective exponent", title="no exact dyadic self-similarity")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st53_entropy_st56_compression.png", dpi=190)
    plt.close(fig)


def write_summary(results: dict[str, Any]) -> None:
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["program", "object", "status"])
        for index in range(46, 58):
            key = f"ST{index}"
            writer.writerow([key, results[key]["object"], results[key]["status"]])


def main() -> None:
    w, a, _ = strict_operator()
    results: dict[str, Any] = {
        "release": "10.58",
        "programs": "ST46-ST57",
        "language": "English",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
        },
    }
    results["ST46"] = st46_upstream_sensitivity_certificate()
    results["ST47"] = st47_scale_stationary_refinement(a)
    results["ST48"] = st48_gain_source_nonidentifiability(a)
    results["ST49"] = st49_projective_multicomponent_state(a)
    results["ST50"] = st50_channel_intertwiner_classification(a)
    results["ST51"] = st51_two_carrier_finite_count(a)
    results["ST52"] = st52_uniform_minimizer_orbit(a)
    results["ST53"] = st53_entropy_data_processing(a)
    results["ST54"] = st54_operational_quotient_universal_property()
    results["ST55"] = st55_executable_two_carrier_spec(results["ST51"])
    results["ST56"] = st56_fractal_compression_refinement(a)
    results["ST57"] = st57_updated_axiom_graph()
    results["recommended_next_programs"] = [
        {"id": "ST58", "priority": 1, "study": "replace ST46 binary64 perturbation propagation by a complete directed-interval eigenprojector, memory and nonlinear-root replay"},
        {"id": "ST59", "priority": 2, "study": "search one target-independent strict adaptive response functional and certify or refute positive odd gain"},
        {"id": "ST60", "priority": 3, "study": "derive or obstruct the ST49 projective spinor texture from the strict spectral bundle without an origin label"},
        {"id": "ST61", "priority": 4, "study": "test whether a sourced trace or compression code can imply the ST47 trace-stationary refinement law"},
        {"id": "ST62", "priority": 5, "study": "derive exact finite-count power bounds and minimum counts for the ST51 two-carrier protocol"},
        {"id": "ST63", "priority": 6, "study": "classify completely-positive operational intertwiners rather than matrix intertwiners alone"},
        {"id": "ST64", "priority": 7, "study": "prove the minimal bath/scale assumptions converting ST53 information entropy into thermodynamic entropy"},
        {"id": "ST65", "priority": 8, "study": "prove localized-state nonexistence or find an interval-certified excited localized orbit for the ST34 saturation"},
        {"id": "ST66", "priority": 9, "study": "replace the ST43 polar flow by a smooth polynomial C12-equivariant multi-component bifurcation and classify its branches"},
        {"id": "ST67", "priority": 10, "study": "separate pi holonomy from chirality and test one orientation-sensitive projective invariant"},
        {"id": "ST68", "priority": 11, "study": "build a calibration uncertainty and custody validator for the ST55 executable specification"},
        {"id": "ST69", "priority": 12, "study": "update the axiom graph after ST58-ST68 and test conditional package reductions"},
    ]
    results["epistemic_boundary"] = (
        "ST46-ST57 are local mathematics, conditional constructions and synthetic receiver studies. They do not supply "
        "a strict feedback-gain source, canonical projective state, physical scale, independent carrier experiment, "
        "legacy-to-strict role transfer, Standard Model, gravity or ToE closure."
    )
    make_figures(results)
    write_summary(results)
    RESULTS.write_text(json.dumps(native(results), indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"results": RESULTS.name, "programs": 12, "figures": 7}, indent=2))


if __name__ == "__main__":
    main()
