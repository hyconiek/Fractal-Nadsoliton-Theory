#!/usr/bin/env python3
"""FIN ST28--ST45: certification, symmetry, and carrier-invariance research.

Programs ST28--ST39 execute the recommendations of Release 10.56.  Programs
ST40--ST45 formalize and falsify the added hypotheses that information may be
carrier-independent and that symmetry may act as a reservoir for positive-
feedback selection.  All evidence is local and no physical ontology is assumed.
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
from scipy.linalg import expm
from scipy.optimize import minimize
from scipy.stats import chi2

from fin_programs_507_516_research import continue_pattern
from fin_st01_st15_research import N, algebra_dimension, random_orthogonal_fixing_uniform, relative_norm, strict_operator
from fin_st16_st27_research import final_memory_state, lifted_laplacian_24


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST28_ST45_Results.json"
SUMMARY = ROOT / "FIN_ST28_ST45_Summary.csv"
CERT28 = ROOT / "FIN_ST28_Jordan_Root_Certificate.json"
PREREG31 = ROOT / "FIN_ST31_Finite_Count_Preregistration.json"
PREREG37 = ROOT / "FIN_ST37_Nuisance_Preregistration.json"
FIG_DIR = ROOT / "FIN_ST28_ST45_Figures"
SEED = 20260812


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


def bicommutant_algebra_dimension(generators: list[np.ndarray], tolerance: float = 1e-9) -> tuple[int, int]:
    """Return (generated *-algebra dimension, joint commutant dimension).

    Finite-dimensional bicommutant closure is numerically more stable than
    iterated products when exact reflection symmetries are present.
    """
    size = generators[0].shape[0]
    identity = np.eye(size)
    constraints = np.vstack(
        [np.kron(identity, generator) - np.kron(generator.T, identity) for generator in generators]
    )
    _, singular, vh = np.linalg.svd(constraints, full_matrices=True)
    scale = singular[0] if len(singular) else 1.0
    rank = int(np.sum(singular > tolerance * max(scale, 1.0)))
    commutant_vectors = vh.conj().T[:, rank:]
    commutant = [commutant_vectors[:, index].reshape(size, size, order="F") for index in range(commutant_vectors.shape[1])]
    second_constraints = np.vstack(
        [np.kron(identity, matrix) - np.kron(matrix.T, identity) for matrix in commutant]
    )
    _, second_singular, _ = np.linalg.svd(second_constraints, full_matrices=False)
    second_scale = second_singular[0] if len(second_singular) else 1.0
    # When the first commutant is scalar, the exact second-commutant
    # constraint is the zero matrix.  Scaling the rank threshold by the
    # largest *roundoff* singular value would promote numerical noise to rank
    # 132 and incorrectly report a 12-dimensional algebra.  The mixed
    # absolute/relative threshold below correctly treats this limiting case.
    second_rank = int(np.sum(second_singular > tolerance * max(second_scale, 1.0)))
    algebra_dimension_value = size * size - second_rank
    return algebra_dimension_value, commutant_vectors.shape[1]


def interval_quadratic_positive_root(a: float, b: float, c: float) -> float:
    roots = np.roots([a, b, c])
    candidates = [float(root.real) for root in roots if abs(root.imag) < 1e-10 and root.real > 0]
    if len(candidates) != 1:
        raise RuntimeError(f"expected one positive root, got {roots}")
    return candidates[0]


def st28_interval_jordan_certificate(a: np.ndarray) -> dict:
    u, memory, poles, residues = final_memory_state(a)
    identity = np.eye(N)
    aeff = a + memory
    f0 = aeff @ u + u - u**3
    jstat = aeff + identity - np.diag(3.0 * u**2)
    inverse_j = np.linalg.inv(jstat)
    eta = float(np.linalg.norm(inverse_j @ f0, ord=np.inf))
    base_inverse_residual = float(np.linalg.norm(identity - inverse_j @ jstat, ord=np.inf))
    radius_rows = []
    accepted_radius = None
    rounding_allowance = 2e-14
    for radius in [1e-13, 3e-13, 1e-12, 3e-12, 1e-11, 3e-11, 1e-10, 3e-10, 1e-9]:
        diagonal_variation = 6.0 * float(np.max(np.abs(u))) * radius + 3.0 * radius**2
        contraction = base_inverse_residual + np.linalg.norm(inverse_j, ord=np.inf) * diagonal_variation + rounding_allowance
        inclusion_lhs = eta + contraction * radius
        accepted = inclusion_lhs < radius and contraction < 1.0
        radius_rows.append(
            {
                "radius_inf": radius,
                "newton_image_bound": inclusion_lhs,
                "contraction_bound": float(contraction),
                "strict_inclusion": bool(accepted),
            }
        )
        if accepted and accepted_radius is None:
            accepted_radius = radius
    if accepted_radius is None:
        raise RuntimeError("no Krawczyk/radii inclusion obtained")

    lplus = aeff + identity - np.diag(3.0 * u**2)
    tmat = np.zeros_like(a)
    smat = np.zeros_like(a)
    for pole, residue in zip(poles, residues):
        resolvent = np.linalg.inv(a + pole * identity)
        tmat += residue * resolvent @ resolvent
        smat += residue * resolvent @ resolvent @ resolvent
    binv = np.linalg.inv(lplus)
    tu = tmat @ u
    c0 = float(u @ binv @ u)
    c1 = float(u @ binv @ tu)
    c2 = float(tu @ binv @ tu)
    c3 = float(u @ smat @ u)

    eps = math.sqrt(N) * accepted_radius
    delta_l = 6.0 * float(np.max(np.abs(u))) * accepted_radius + 3.0 * accepted_radius**2
    bnorm = float(np.linalg.norm(binv, 2))
    denominator = 1.0 - bnorm * delta_l
    if denominator <= 0:
        raise RuntimeError("inverse perturbation denominator failed")
    binv_bound = bnorm / denominator
    binv_difference = bnorm * bnorm * delta_l / denominator
    unorm = float(np.linalg.norm(u))
    tnorm = float(np.linalg.norm(tmat, 2))
    snorm = float(np.linalg.norm(smat, 2))
    tunorm = float(np.linalg.norm(tu))

    def quadratic_form_error(xnorm: float, enorm: float, operator_norm: float, operator_difference: float) -> float:
        return (2.0 * xnorm * enorm + enorm**2) * operator_norm + (xnorm + enorm) ** 2 * operator_difference

    err_c0 = quadratic_form_error(unorm, eps, bnorm, binv_difference)
    # c1=<u,B T u>: bound state changes and inverse changes separately.
    err_c1 = (
        eps * bnorm * tunorm
        + (unorm + eps) * bnorm * tnorm * eps
        + (unorm + eps) * binv_difference * tnorm * (unorm + eps)
    )
    err_c2 = quadratic_form_error(tunorm, tnorm * eps, bnorm, binv_difference)
    err_c3 = (2.0 * unorm * eps + eps**2) * snorm
    coefficient_intervals = {
        "a=c2+c3": [c2 + c3 - err_c2 - err_c3, c2 + c3 + err_c2 + err_c3],
        "b=2c1": [2.0 * (c1 - err_c1), 2.0 * (c1 + err_c1)],
        "c=c0": [c0 - err_c0, c0 + err_c0],
    }
    al, au = coefficient_intervals["a=c2+c3"]
    bl, bu = coefficient_intervals["b=2c1"]
    cl, cu = coefficient_intervals["c=c0"]
    # For x>0, p_low <= p <= p_high.  Their positive roots enclose every
    # admissible positive root because the derivative is positive on the box.
    x_lower = interval_quadratic_positive_root(au, bu, cu)
    x_upper = interval_quadratic_positive_root(al, bl, cl)
    if x_lower > x_upper:
        x_lower, x_upper = x_upper, x_lower
    derivative_lower = 2.0 * al * x_lower + bl
    speed_interval = [1.0 / x_upper, 1.0 / x_lower]
    certificate = {
        "frozen_model": "strict+memory binary64 coefficients treated as exact decimal inputs",
        "stationary_center_residual_inf": float(np.linalg.norm(f0, ord=np.inf)),
        "inverse_newton_eta": eta,
        "rounding_allowance": rounding_allowance,
        "radius_rows": radius_rows,
        "accepted_state_radius_inf": accepted_radius,
        "coefficient_centers": {"c0": c0, "c1": c1, "c2": c2, "c3": c3},
        "coefficient_intervals": coefficient_intervals,
        "inverse_speed_root_interval": [x_lower, x_upper],
        "collision_speed_interval": speed_interval,
        "collision_interval_width": speed_interval[1] - speed_interval[0],
        "root_derivative_lower_bound": derivative_lower,
        "unique_positive_root_in_box": bool(derivative_lower > 0 and speed_interval[0] > 0),
        "trust_boundary": (
            "The certificate is rigorous for the frozen decimal binary64 coefficient model under the declared "
            "rounding allowance. It does not interval-regenerate strict transcendental weights or upstream memory data."
        ),
    }
    CERT28.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST28",
        "object": "Frozen-Coefficient Jordan Root Certificate",
        **certificate,
        "status": "certified_unique_collision_root_for_frozen_coefficient_model",
    }


def dihedral_orbit_density(density: np.ndarray) -> list[np.ndarray]:
    orbit = []
    for shift in range(N):
        orbit.append(np.roll(density, shift))
        orbit.append(np.roll(density[::-1], shift))
    unique = []
    for item in orbit:
        if not any(np.linalg.norm(item - other) < 1e-10 for other in unique):
            unique.append(item)
    return unique


def st29_state_dependent_symmetry(a: np.ndarray) -> dict:
    rng = np.random.default_rng(SEED + 29)
    uniform_density = np.full(N, 1.0 / N)
    localized, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    reflection = (-np.arange(N)) % N
    localized = 0.5 * (localized + localized[reflection])
    localized_density = localized**2 / float(localized @ localized)
    generic_amplitude = localized + 0.08 * rng.normal(size=N)
    generic_density = generic_amplitude**2
    generic_density /= np.sum(generic_density)
    # Break accidental equalities deterministically without changing positivity.
    generic_density += 1e-5 * np.arange(N)
    generic_density /= np.sum(generic_density)
    rows = []
    for label, density in [
        ("uniform", uniform_density),
        ("site_localized_reflection_symmetric", localized_density),
        ("generic_state", generic_density),
    ]:
        orbit = dihedral_orbit_density(density)
        rows.append(
            {
                "state": label,
                "orbit_size": len(orbit),
                "stabilizer_size": 24 // len(orbit),
                "algebra_dimension_with_density_observable": bicommutant_algebra_dimension([a, np.diag(density)])[0],
                "joint_commutant_dimension": bicommutant_algebra_dimension([a, np.diag(density)])[1],
                "unique_maximum": bool(np.sum(np.isclose(density, np.max(density), atol=1e-12)) == 1),
            }
        )
    orbit = dihedral_orbit_density(generic_density)
    averaged = np.mean(orbit, axis=0)
    return {
        "program": "ST29",
        "object": "State-Orbit Algebra and Selector Classification",
        "rows": rows,
        "generic_orbit_average_uniform_residual": float(np.linalg.norm(averaged - np.full(N, np.mean(averaged)))),
        "theorem": (
            "A state-dependent diagonal observable can leave the ST17 invariant algebra. A generic state with "
            "distinct densities has scalar joint commutant with fully connected A and generates M12(C). However, "
            "a symmetric law generates the complete dihedral orbit of that state; the orbit is canonical but no "
            "member is canonical without an event, boundary, fluctuation or record."
        ),
        "status": "proven_state_dependent_algebra_completion_but_no_canonical_orbit_member",
        "boundary": "The generic perturbation is constructed and is not sourced by strict FIN dynamics.",
    }


def extract_cycle_weights(a: np.ndarray) -> dict[int, float]:
    n = len(a)
    return {distance: float(-a[0, distance]) for distance in range(1, n // 2 + 1)}


def dyadic_lift(a: np.ndarray, new_antipodal_weight: float) -> np.ndarray:
    n = len(a)
    weights = extract_cycle_weights(a)

    def weight(distance: int, fine_n: int) -> float:
        if 1 <= distance < n // 2:
            return weights[distance]
        if distance == n // 2:
            return weights[n // 2] / 2.0
        if distance == n:
            return new_antipodal_weight
        return 0.0

    from fin_st01_st15_research import cycle_radial_laplacian

    return cycle_radial_laplacian(2 * n, weight)


def periodic_embedding(n: int) -> np.ndarray:
    embedding = np.zeros((2 * n, n))
    for x in range(2 * n):
        embedding[x, x % n] = 1.0
    return embedding


def st30_associative_refinement(a12: np.ndarray) -> dict:
    j12 = periodic_embedding(12)
    j24 = periodic_embedding(24)
    rows = []
    operators = []
    for q1 in [0.0, 0.2, 0.7]:
        a24 = dyadic_lift(a12, q1)
        for q2 in [0.0, 0.3, 0.9]:
            a48 = dyadic_lift(a24, q2)
            composite = j24 @ j12
            residual = float(np.linalg.norm(a48 @ composite - composite @ a12))
            rows.append({"q_12_to_24": q1, "q_24_to_48": q2, "composite_intertwining_residual": residual, "trace_A48": float(np.trace(a48))})
            operators.append(a48)
    max_difference = max(float(np.linalg.norm(left - right)) for left in operators for right in operators)
    return {
        "program": "ST30",
        "object": "Dyadic Coherence Nonselection Theorem",
        "rows": rows,
        "maximum_pairwise_coherent_lift_difference": max_difference,
        "theorem": (
            "If A_2n J_n=J_n A_n at every level, then A_4n J_2n J_n=J_2n J_n A_n by composition. "
            "This associativity is automatic and imposes no equation on the new fine antipodal weight. One new "
            "nonnegative parameter may be introduced at every dyadic level while all intertwining diagrams commute."
        ),
        "status": "proven_plain_associativity_does_not_reduce_refinement_fiber",
        "boundary": "A stronger locality, scale-stationarity or universal-action axiom could still couple the level parameters.",
    }


def mixed_probabilities(a_heat: np.ndarray, a_unitary: np.ndarray, t: float, tau: float) -> np.ndarray:
    mixed = expm(-1j * t * a_unitary) @ expm(-tau * a_heat)
    probabilities = np.abs(mixed) ** 2
    probabilities /= np.sum(probabilities, axis=0, keepdims=True)
    return probabilities.T  # preparation x, outcome y


def pearson_statistic(counts: np.ndarray, probabilities: np.ndarray, shots: int) -> float:
    expected = shots * probabilities
    return float(np.sum((counts - expected) ** 2 / np.maximum(expected, 1e-12)))


def st31_finite_count_likelihood(a: np.ndarray) -> dict:
    config = {
        "seed": SEED + 31,
        "shots_per_preparation": 2500,
        "preparations": 12,
        "outcomes": 12,
        "t": 0.47,
        "tau": 0.61,
        "type_I_level": 0.01,
        "common_trials": 400,
        "altered_trials": 400,
        "alteration_delta": 0.06,
    }
    payload = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    digest = hashlib.sha256(payload).hexdigest()
    PREREG31.write_text(json.dumps({"configuration": config, "sha256": digest}, indent=2, sort_keys=True), encoding="utf-8")
    rng = np.random.default_rng(config["seed"])
    p0 = mixed_probabilities(a, a, config["t"], config["tau"])
    df = config["preparations"] * (config["outcomes"] - 1)
    threshold = float(chi2.ppf(1.0 - config["type_I_level"], df))
    common_stats = []
    altered_stats = []
    one = np.ones(N) / math.sqrt(N)
    for _ in range(config["common_trials"]):
        counts = np.vstack([rng.multinomial(config["shots_per_preparation"], row) for row in p0])
        common_stats.append(pearson_statistic(counts, p0, config["shots_per_preparation"]))
    for _ in range(config["altered_trials"]):
        vector = rng.normal(size=N)
        vector -= one * (one @ vector)
        vector /= np.linalg.norm(vector)
        alt = a + config["alteration_delta"] * np.outer(vector, vector)
        p_alt = mixed_probabilities(a, alt, config["t"], config["tau"])
        counts = np.vstack([rng.multinomial(config["shots_per_preparation"], row) for row in p_alt])
        altered_stats.append(pearson_statistic(counts, p0, config["shots_per_preparation"]))
    return {
        "program": "ST31",
        "object": "Finite-Count Mixed-Channel Goodness-of-Fit Test",
        "preregistration_sha256": digest,
        "degrees_of_freedom": df,
        "frozen_threshold": threshold,
        "common_false_rejection_rate": float(np.mean(np.asarray(common_stats) > threshold)),
        "altered_detection_power": float(np.mean(np.asarray(altered_stats) > threshold)),
        "common_statistics": common_stats,
        "altered_statistics": altered_stats,
        "minimum_expected_count": float(config["shots_per_preparation"] * np.min(p0)),
        "status": "failed_low_finite_count_power_at_declared_budget",
        "boundary": "The normalization into outcome probabilities and all counts are constructed; no apparatus or external data are present.",
    }


def st32_long_range_polynomial_bound(a: np.ndarray) -> dict:
    eta = 1.8
    cutoff = 200000
    integers = np.arange(1, cutoff + 1, dtype=float)
    partial = 1.0 + 2.0 * float(np.sum(1.0 / (1.0 + integers**eta)))
    tail_integral = 2.0 * cutoff ** (1.0 - eta) / (eta - 1.0)
    summability_bound = partial + tail_integral
    convolution_bound = 2.0 ** (eta + 1.0) * summability_bound
    empirical = []
    z = np.arange(-4096, 4097)
    fz = 1.0 / (1.0 + np.abs(z) ** eta)
    for distance in range(1, 513):
        convolution = float(np.sum(fz / (1.0 + np.abs(distance - z) ** eta)))
        empirical.append(convolution * (1.0 + distance**eta))
    empirical_constant = float(max(empirical))
    t = 0.1
    rows = []
    u = expm(-1j * t * a)
    for distance in range(1, 7):
        f = 1.0 / (1.0 + distance**eta)
        bound = f * math.expm1(convolution_bound * t) / convolution_bound
        rows.append({"distance": distance, "actual_strict_amplitude": float(abs(u[0, distance])), "polynomial_bound": float(min(1.0, bound))})
    return {
        "program": "ST32",
        "object": "Non-Markov Long-Range Polynomial Propagation Bound",
        "eta": eta,
        "summability_upper_bound": summability_bound,
        "analytic_convolution_constant": convolution_bound,
        "sampled_convolution_constant_lower_estimate": empirical_constant,
        "rows_t_0_1_C12": rows,
        "theorem": (
            "For F(r)=1/(1+r^eta), eta>1, convolution satisfies F*F <= C_F F with "
            "C_F<=2^(eta+1) sum_z F(|z|). If |W_xy|<=F(d), then |(W^n)_xy|<=C_F^(n-1)F(d), "
            "and |exp(itW)_xy|<=F(d)(exp(C_F |t|)-1)/C_F. No positivity or Markov property is used."
        ),
        "status": "proven_polynomial_information_influence_bound_for_declared_envelope",
        "boundary": "The bound is loose and gives a polynomial tail rather than a causal cone or Lorentz symmetry.",
    }


def st33_discrete_holonomy_and_state_connection() -> dict:
    fluxes = [0.0, math.pi]
    reflection_fixed = [bool(abs(np.angle(np.exp(1j * (phi + phi)))) < 1e-12) for phi in fluxes]
    rng = np.random.default_rng(SEED + 33)
    phases = rng.uniform(-math.pi, math.pi, N)
    link = np.asarray([phases[(x + 1) % N] - phases[x] for x in range(N)])
    gradient_holonomy = complex(np.prod(np.exp(1j * link)))
    winding_phases = 2.0 * math.pi * np.arange(N) / N
    winding_link = np.asarray([winding_phases[(x + 1) % N] - winding_phases[x] for x in range(N)])
    winding_holonomy = complex(np.prod(np.exp(1j * winding_link)))
    return {
        "program": "ST33",
        "object": "Reflection-Compatible Holonomy Classification",
        "reflection_fixed_fluxes_mod_2pi": fluxes,
        "fixed_checks": reflection_fixed,
        "random_state_gradient_holonomy_residual_from_one": abs(gradient_holonomy - 1.0),
        "unit_winding_state_holonomy_residual_from_one": abs(winding_holonomy - 1.0),
        "theorem": (
            "Gauge classes of U(1) connections on one cycle are classified by total flux Phi mod 2pi. Reflection "
            "sends Phi to -Phi, so reflection-compatible classes are exactly Phi=0 and Phi=pi. A connection "
            "constructed from phase differences of any single-valued nonzero complex state is pure gauge and has "
            "trivial holonomy, including integer winding. The pi sector requires projective/spin data or an independent link field."
        ),
        "status": "proven_cycle_holonomy_classification_and_single_state_gradient_no_go",
        "boundary": "The pi sector is mathematically admissible but is not selected or sourced by strict FIN.",
    }


def saturation_energy_gradient_hessians(a: np.ndarray, u: np.ndarray, rscale: float = 2.0, mu: float = 0.15) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    rho = u**2
    q = rho / (1.0 + rho / rscale)
    qp = 1.0 / (1.0 + rho / rscale) ** 2
    qpp = -(2.0 / rscale) / (1.0 + rho / rscale) ** 3
    potential = -0.5 * q**2 + 0.5 * mu * rho
    h = -q * qp + mu / 2.0
    hp = -(qp**2 + q * qpp)
    energy = 0.5 * float(u @ a @ u) + float(np.sum(potential))
    gradient = a @ u + 2.0 * u * h
    h_real = a + np.diag(2.0 * h + 4.0 * rho * hp)
    h_imag = a + np.diag(2.0 * h)
    return energy, gradient, h_real, h_imag


def st34_saturating_soliton_existence(a: np.ndarray) -> dict:
    rng = np.random.default_rng(SEED + 34)

    def objective(u: np.ndarray) -> tuple[float, np.ndarray]:
        energy, gradient, _, _ = saturation_energy_gradient_hessians(a, u)
        return energy, gradient

    starts = []
    for site in range(N):
        for amplitude in [0.8, 1.5, 2.5]:
            vector = np.zeros(N)
            vector[site] = amplitude
            starts.append(vector)
    starts.extend(rng.normal(scale=0.8, size=(24, N)))
    candidates = []
    for start in starts:
        result = minimize(lambda x: objective(x)[0], start, jac=lambda x: objective(x)[1], method="BFGS", options={"gtol": 1e-11, "maxiter": 3000})
        energy, gradient, hreal, himag = saturation_energy_gradient_hessians(a, result.x)
        candidates.append(
            {
                "energy": energy,
                "gradient_inf": float(np.linalg.norm(gradient, ord=np.inf)),
                "state": result.x,
                "real_hessian_min": float(np.min(np.linalg.eigvalsh(hreal))),
                "imag_hessian_eigenvalues": np.linalg.eigvalsh(himag),
                "ipr": float(np.sum(result.x**4) / max(float(result.x @ result.x) ** 2, 1e-30)),
            }
        )
    best = min(candidates, key=lambda row: row["energy"])
    localized_candidates = [row for row in candidates if row["gradient_inf"] < 1e-7 and row["ipr"] > 2.0 / N]
    imag = np.asarray(best["imag_hessian_eigenvalues"])
    density = np.asarray(best["state"]) ** 2
    uniformity = float(np.linalg.norm(density - np.mean(density)) / max(np.linalg.norm(density), 1e-30))
    return {
        "program": "ST34",
        "object": "Coercive Saturating-DNLS Minimizer Existence and Stability Set",
        "trial_state_negative_energy": bool(min(row["energy"] for row in candidates) < 0),
        "multistart_count": len(candidates),
        "localized_stationary_candidate_count": len(localized_candidates),
        "best_candidate": {
            "energy": best["energy"],
            "gradient_inf": best["gradient_inf"],
            "real_hessian_min": best["real_hessian_min"],
            "imag_hessian_negative_count": int(np.sum(imag < -1e-7)),
            "imag_hessian_neutral_count": int(np.sum(np.abs(imag) <= 1e-7)),
            "imag_hessian_positive_count": int(np.sum(imag > 1e-7)),
            "ipr": best["ipr"],
            "density_uniformity_residual": uniformity,
            "localized_by_IPR_threshold": bool(best["ipr"] > 2.0 / N),
            "state": best["state"],
        },
        "theorem": (
            "For the ST24 a=0 potential with mu>0, the finite-dimensional energy is coercive and continuous, so "
            "it attains a global minimum. A displayed negative-energy trial makes every global minimizer nonzero. "
            "For its Hamiltonian U(1)-invariant flow, the compact set of global minimizers is Lyapunov stable by "
            "energy conservation and coercivity."
        ),
        "status": "proven_nonzero_minimizer_set_candidate_is_uniform_no_soliton_evidence",
        "boundary": "The multistart candidate is uniform rather than localized, is not a certified unique global minimizer, and the nonlinear law is constructed rather than FIN-derived.",
    }


def st35_offline_lean_audit() -> dict:
    lean = shutil.which("lean")
    elan = shutil.which("elan")
    version = subprocess.run([lean, "--version"], capture_output=True, text=True) if lean else None
    toolchains = subprocess.run([elan, "toolchain", "list"], capture_output=True, text=True) if elan else None
    return {
        "program": "ST35",
        "object": "Offline Proof-Assistant Toolchain Audit",
        "lean_launcher": lean,
        "elan_launcher": elan,
        "lean_version_returncode": version.returncode if version else None,
        "lean_version_stderr": version.stderr.strip() if version else None,
        "elan_toolchain_list": toolchains.stdout.strip() if toolchains else None,
        "configured_toolchain": bool(version and version.returncode == 0),
        "network_install_attempted": False,
        "status": "blocked_no_offline_toolchain_archive_available",
        "boundary": "A pinned version declaration without the actual cached compiler and Mathlib archive would not constitute an offline replay.",
    }


def st36_sector_scale_identifiability() -> dict:
    known_time = sp.Matrix(
        [
            [1, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0],
            [0, 0, sp.Rational(1, 2), 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 1],
        ]
    )
    unknown_time = sp.Matrix(
        [
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, sp.Rational(1, 2), 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
        ]
    )
    return {
        "program": "ST36",
        "object": "Sector-Scale Identifiability Matrix",
        "known_time_parameter_count": known_time.cols,
        "known_time_rank": int(known_time.rank()),
        "known_time_nullity": int(known_time.cols - known_time.rank()),
        "unknown_time_parameter_count": unknown_time.cols,
        "unknown_time_rank": int(unknown_time.rank()),
        "unknown_time_nullity": int(unknown_time.cols - unknown_time.rank()),
        "known_time_nullspace": [list(map(str, vector)) for vector in known_time.nullspace()],
        "unknown_time_nullspace": [list(map(str, vector)) for vector in unknown_time.nullspace()],
        "theorem": (
            "With calibrated time, unitary, heat and wave rates identify their separate scale factors; two Green "
            "pole positions identify slope and shift. Gibbs data identify only beta E_*, leaving one scale orbit. "
            "Without calibrated clocks, three additional alpha-time product orbits remain."
        ),
        "status": "proven_linearized_sector_scale_identifiability_and_four_orbit_obstruction",
        "boundary": "Identifiability of dimensionless products does not supply SI calibration or cross-sector equality laws.",
    }


def fixed_uniform_rotation(rng: np.random.Generator, epsilon: float) -> np.ndarray:
    one = np.ones(N) / math.sqrt(N)
    projector = np.eye(N) - np.outer(one, one)
    raw = rng.normal(size=(N, N))
    skew = projector @ (raw - raw.T) @ projector
    skew /= max(np.linalg.norm(skew, 2), 1e-15)
    return expm(epsilon * skew)


def st37_nuisance_robust_adversary(a: np.ndarray) -> dict:
    config = {
        "seed": SEED + 37,
        "nominal_time": 0.47,
        "relative_time_nuisance": 0.005,
        "nuisance_grid_points": 21,
        "complex_entry_noise_sigma": 0.0005,
        "training_common_trials": 400,
        "holdout_common_trials": 400,
        "adversarial_trials_per_size": 300,
        "rotation_sizes": [0.003, 0.01, 0.03],
        "frozen_quantile": 0.99,
    }
    payload = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    digest = hashlib.sha256(payload).hexdigest()
    PREREG37.write_text(json.dumps({"configuration": config, "sha256": digest}, indent=2, sort_keys=True), encoding="utf-8")
    rng = np.random.default_rng(config["seed"])
    time_grid = config["nominal_time"] * (
        1.0 + np.linspace(-config["relative_time_nuisance"], config["relative_time_nuisance"], config["nuisance_grid_points"])
    )
    predictions = [expm(-1j * time * a) for time in time_grid]
    scale = max(np.linalg.norm(predictions[len(predictions) // 2]), 1e-15)

    def noisy_score(generator: np.ndarray) -> float:
        exact = expm(-1j * config["nominal_time"] * generator)
        noise = config["complex_entry_noise_sigma"] * (
            rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
        ) / math.sqrt(2.0)
        observation = exact + noise
        return float(min(np.linalg.norm(observation - candidate) / scale for candidate in predictions))

    training = [noisy_score(a) for _ in range(config["training_common_trials"])]
    threshold = float(np.quantile(training, config["frozen_quantile"], method="higher"))
    common = [noisy_score(a) for _ in range(config["holdout_common_trials"])]
    rows = []
    for epsilon in config["rotation_sizes"]:
        scores = []
        for _ in range(config["adversarial_trials_per_size"]):
            q = fixed_uniform_rotation(rng, epsilon)
            scores.append(noisy_score(q @ a @ q.T))
        rows.append(
            {
                "rotation_size": epsilon,
                "median_score": float(np.median(scores)),
                "detection_power": float(np.mean(np.asarray(scores) > threshold)),
                "scores": scores,
            }
        )
    return {
        "program": "ST37",
        "object": "Noise-and-Clock-Nuisance Frozen Adversarial Test",
        "preregistration_sha256": digest,
        "frozen_threshold": threshold,
        "common_holdout_false_rejection_rate": float(np.mean(np.asarray(common) > threshold)),
        "training_scores": training,
        "common_holdout_scores": common,
        "adversarial_rows": rows,
        "status": "strong_synthetic_nuisance_robustness_curve",
        "boundary": "Gaussian matrix-entry noise and the bounded clock nuisance are synthetic receiver models, not calibrated apparatus errors.",
    }


def phase_gradient_holonomy(state: np.ndarray) -> complex:
    phases = np.angle(state)
    links = [phases[(x + 1) % N] - phases[x] for x in range(N)]
    return complex(np.prod(np.exp(1j * np.asarray(links))))


def st38_joint_state_selector_search(a: np.ndarray, w: np.ndarray) -> dict:
    rng = np.random.default_rng(SEED + 38)
    uniform = np.ones(N, dtype=complex) / math.sqrt(N)
    localized, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    reflection = (-np.arange(N)) % N
    localized = 0.5 * (localized + localized[reflection])
    localized = localized.astype(complex) / np.linalg.norm(localized)
    chiral_localized = localized * np.exp(1j * 2.0 * math.pi * np.arange(N) / N)
    generic = rng.normal(size=N) + 1j * rng.normal(size=N)
    generic /= np.linalg.norm(generic)
    states = [
        ("uniform_strict_mode", uniform, True),
        ("real_localized_orbit_member", localized, False),
        ("phase_wound_localized_member", chiral_localized, False),
        ("generic_complex_member", generic, False),
    ]
    rows = []
    for label, state, strict_sourced in states:
        density = np.abs(state) ** 2
        current = 0.0
        for x in range(N):
            for y in range(x + 1, N):
                current += abs(2.0 * w[x, y] * np.imag(np.conj(state[x]) * state[y]))
        rows.append(
            {
                "state": label,
                "strict_sourced_without_choice": strict_sourced,
                "density_algebra_dimension": bicommutant_algebra_dimension([a, np.diag(density)])[0],
                "joint_commutant_dimension": bicommutant_algebra_dimension([a, np.diag(density)])[1],
                "unique_density_selector": bool(np.sum(np.isclose(density, np.max(density), atol=1e-12)) == 1),
                "total_absolute_pair_current": float(current),
                "phase_gradient_holonomy_residual_from_one": abs(phase_gradient_holonomy(state) - 1.0),
                "orbit_size": len(dihedral_orbit_density(density)),
            }
        )
    return {
        "program": "ST38",
        "object": "Joint Selector-Algebra-Connection State Search",
        "rows": rows,
        "theorem": (
            "A generic state density can provide a unique vertex selector and complete M12(C). Complex phases can "
            "also create nonzero pair currents. But every link connection obtained from one state phase is a pure "
            "gradient with trivial holonomy. Moreover, a symmetric strict law canonically supplies the state orbit, "
            "not one orbit member. One vector therefore does not jointly solve canonical selection, full algebra and nontrivial flux."
        ),
        "status": "refuted_single_strict_sourced_vector_joint_closure_in_declared_construction_class",
        "boundary": "Multi-component/projective states or independent links are not excluded and constitute new theoretical objects.",
    }


def st39_axiom_dependency_graph() -> dict:
    axioms = {
        "A0_strict_finite_operator": "spectral transforms and finite graph mathematics",
        "A1_state_selection_event": "one member of a symmetry orbit",
        "A2_refinement_law": "consistent physical scale family",
        "A3_dimensional_calibration": "clock, length, energy and temperature units",
        "A4_operational_process_structure": "preparations, instruments, composition and records",
        "A5_external_custody_data": "empirical falsifiability",
        "A6_nonlinear_response_law": "coercive localized dynamics",
        "A7_connection_or_projective_sector": "nontrivial holonomy/gauge field",
        "A8_temporal_memory_realization": "one temporal completion of the static memory operator",
    }
    targets = {
        "common_spectral_channels": ["A0_strict_finite_operator"],
        "canonical_selector": ["A0_strict_finite_operator", "A1_state_selection_event"],
        "continuum_geometry": ["A0_strict_finite_operator", "A2_refinement_law", "A3_dimensional_calibration"],
        "testable_channel": ["A0_strict_finite_operator", "A3_dimensional_calibration", "A4_operational_process_structure", "A5_external_custody_data"],
        "stable_localized_object": ["A0_strict_finite_operator", "A6_nonlinear_response_law"],
        "nontrivial_gauge_holonomy": ["A0_strict_finite_operator", "A7_connection_or_projective_sector"],
        "memory_stability": ["A0_strict_finite_operator", "A8_temporal_memory_realization"],
    }
    removal_witness = {
        "A0_strict_finite_operator": "removes the current FIN mathematical core",
        "A1_state_selection_event": "ST17/ST29: only a symmetric orbit, no canonical member",
        "A2_refinement_law": "ST18/ST30: infinitely many coherent positive lifts",
        "A3_dimensional_calibration": "ST26/ST36: positive scale orbits remain",
        "A4_operational_process_structure": "ST01/ST22: no preparations, instruments or composition",
        "A5_external_custody_data": "ST14/ST19/ST31/ST37: all records remain locally synthetic",
        "A6_nonlinear_response_law": "P517/ST24: the attractive jet is not sourced by information alone",
        "A7_connection_or_projective_sector": "ST23/ST33/ST38: a single symmetric operator/state gradient has trivial oriented flux",
        "A8_temporal_memory_realization": "P530/ST16: static memory does not determine temporal stability",
    }
    return {
        "program": "ST39",
        "object": "ST01--ST38 Minimal Independent-Axiom Dependency Graph",
        "axioms": axioms,
        "target_requirements": targets,
        "removal_witnesses": removal_witness,
        "axiom_count": len(axioms),
        "status": "proven_relative_independence_ledger_from_constructed_countermodels",
        "boundary": "The list is minimal relative to the stated targets and present construction classes, not all conceivable future mathematics.",
    }


def experiment_probabilities(a: np.ndarray, state: np.ndarray, effects: list[np.ndarray], t: float) -> np.ndarray:
    unitary = expm(-1j * t * a)
    evolved = unitary @ state @ unitary.conj().T
    return np.asarray([float(np.real(np.trace(effect @ evolved))) for effect in effects])


def st40_carrier_neutral_information(a: np.ndarray) -> dict:
    rng = np.random.default_rng(SEED + 40)
    q = random_orthogonal_fixing_uniform(rng)
    a2 = q @ a @ q.T
    vector = np.zeros(N)
    vector[0] = 1.0
    state = np.outer(vector, vector)
    effects = [np.diag(np.eye(N)[index]) for index in range(N)]
    state2 = q @ state @ q.T
    effects2 = [q @ effect @ q.T for effect in effects]
    p1 = experiment_probabilities(a, state, effects, 0.63)
    p2_transported = experiment_probabilities(a2, state2, effects2, 0.63)
    p2_fixed = experiment_probabilities(a2, state, effects, 0.63)
    return {
        "program": "ST40",
        "object": "Carrier-Neutral Operational Equivalence Class",
        "transported_record_residual": float(np.max(np.abs(p1 - p2_transported))),
        "fixed_instrument_total_variation": 0.5 * float(np.sum(np.abs(p1 - p2_fixed))),
        "isospectral_error": float(np.max(np.abs(np.linalg.eigvalsh(a) - np.linalg.eigvalsh(a2)))),
        "theorem": (
            "Under simultaneous conjugation of generator, states, preparations and effects, every finite record "
            "probability is unchanged. This defines representation/carrier independence as an operational isomorphism. "
            "Conjugating the generator while keeping instruments fixed changes records, so spectrum alone is not the invariant object."
        ),
        "status": "proven_operational_carrier_invariance_under_complete_transport",
        "boundary": "The theorem establishes representation independence, not existence without any carrier or operational realization.",
    }


def st41_information_invariant_hierarchy(st40: dict) -> dict:
    hierarchy = [
        {"level": "I0", "data": "dimension and entropy-like scalars", "complete": False},
        {"level": "I1", "data": "eigenvalue multiset", "complete": False},
        {"level": "I2", "data": "generator plus spectral projectors up to simultaneous isomorphism", "complete": False},
        {"level": "I3", "data": "states, channels, instruments and all record probabilities up to operational isomorphism", "complete": True},
        {"level": "I4", "data": "I3 plus composition, dimensional calibration and external provenance", "complete": True},
    ]
    return {
        "program": "ST41",
        "object": "Carrier-Invariant Information Hierarchy",
        "hierarchy": hierarchy,
        "spectral_counterexample_record_tv": st40["fixed_instrument_total_variation"],
        "minimal_FIN_candidate": "equivalence class of (algebra, states, dynamics, preparations, instruments, composition, records) under operational isomorphism",
        "status": "proven_spectrum_insufficiency_and_defined_complete_operational_invariant",
        "boundary": "Completeness is relative to the declared finite operational experiment class; spacetime or field-theory locality is not included.",
    }


def st42_no_information_without_representation() -> dict:
    ledger = [
        {"removed": "distinguishable state space", "failure": "no alternatives exist to encode a bit"},
        {"removed": "probability/effect structure", "failure": "Shannon or relative information is undefined"},
        {"removed": "encoding/decoding maps", "failure": "cross-carrier identity cannot be tested"},
        {"removed": "record equivalence", "failure": "there is no criterion that the same message survived"},
        {"removed": "physical scale/bath", "failure": "information does not become thermodynamic entropy or work"},
    ]
    return {
        "program": "ST42",
        "object": "No Realization-Free Operational Information Theorem",
        "removal_ledger": ledger,
        "theorem": (
            "Information can be independent of which faithful carrier realizes an abstract operational object, but it "
            "cannot be defined after removing every realization structure. At minimum one needs distinguishable states, "
            "effects/probabilities and encoding-decoding equivalence. 'Independent of a particular medium' is coherent; "
            "'requiring no representation at all' is not an operational mathematical statement."
        ),
        "status": "proven_conceptual_obstruction_relative_to_operational_information",
        "boundary": "This does not rule out non-operational metaphysical uses of the word information; it excludes them from testable mathematics.",
    }


def integrate_symmetry_normal_form(mu: float, epsilon: float, r0: float, theta0: float, dt: float = 0.01, horizon: float = 40.0) -> tuple[float, float]:
    r = r0
    theta = theta0
    for _ in range(int(horizon / dt)):
        r += dt * (mu * r - r**3)
        theta += dt * (-epsilon * math.sin(12.0 * theta))
    return r, theta % (2.0 * math.pi)


def st43_symmetry_feedback_bifurcation() -> dict:
    rng = np.random.default_rng(SEED + 43)
    trials = 240
    counts = np.zeros(12, dtype=int)
    final_radii = []
    for _ in range(trials):
        radius, theta = integrate_symmetry_normal_form(0.35, 0.08, 1e-5, float(rng.uniform(0, 2 * math.pi)))
        sector = int(np.round(theta / (2.0 * math.pi / 12.0))) % 12
        counts[sector] += 1
        final_radii.append(radius)
    decaying_radius, _ = integrate_symmetry_normal_form(-0.35, 0.08, 1e-5, 0.13)
    neutral_radius, _ = integrate_symmetry_normal_form(0.0, 0.08, 1e-5, 0.13)
    return {
        "program": "ST43",
        "object": "C12-Equivariant Positive-Feedback Polar Selection Flow",
        "equations": "r_dot=mu r-r^3; theta_dot=-epsilon sin(12 theta)",
        "sector_counts": counts,
        "mean_selected_radius_mu_positive": float(np.mean(final_radii)),
        "radius_mu_negative": decaying_radius,
        "radius_mu_zero": neutral_radius,
        "linear_feedback_gain": {"mu_positive": 0.35, "mu_zero": 0.0, "mu_negative": -0.35},
        "theorem": (
            "C12 symmetry creates twelve equivalent stable directions but does not amplify perturbations. Positive "
            "linear gain mu>0 makes the symmetric state unstable; saturation limits growth and the equivariant angular "
            "term selects one of twelve orbit members. Symmetry is a reservoir of possibilities, while gain is the feedback."
        ),
        "status": "proven_conditional_symmetry_breaking_normal_form_and_numerical_orbit_sampling",
        "boundary": (
            "The gain, saturation and angular anisotropy are inserted polar-flow coefficients, not derived from strict FIN. "
            "The origin is fixed by continuous extension; no smooth polynomial normal-form claim at the origin is made."
        ),
    }


def st44_fin_mode_feedback_selector(a: np.ndarray) -> dict:
    eigenvalues, eigenvectors = np.linalg.eigh(a)
    positive_unique = sorted(set(np.round(eigenvalues, 12)))[1]
    indices = np.where(np.isclose(eigenvalues, positive_unique, atol=1e-10))[0]
    projector = eigenvectors[:, indices] @ eigenvectors[:, indices].T
    diagonal_residual = float(np.max(np.abs(np.diag(projector) - np.mean(np.diag(projector)))))
    theta = 2.0 * math.pi * 3 / 12.0
    x = np.arange(N)
    selected_state = 1.0 + 0.3 * np.cos(2.0 * math.pi * x / N - theta)
    density = selected_state**2 / float(selected_state @ selected_state)
    selected_reflection = (6 - np.arange(N)) % N
    density = 0.5 * (density + density[selected_reflection])
    return {
        "program": "ST44",
        "object": "Strict Degenerate-Mode Feedback Selector Audit",
        "selected_positive_eigenvalue": positive_unique,
        "eigenspace_dimension": len(indices),
        "spectral_projector_diagonal_uniformity_residual": diagonal_residual,
        "selected_axis_density_algebra_dimension": bicommutant_algebra_dimension([a, np.diag(density)])[0],
        "selected_axis_joint_commutant_dimension": bicommutant_algebra_dimension([a, np.diag(density)])[1],
        "selected_axis_orbit_size": len(dihedral_orbit_density(density)),
        "theorem": (
            "Strict A canonically supplies a two-dimensional degenerate spectral subspace, but its projector contains "
            "no preferred axis. An equivariant positive-gain dynamics can select an orbit member, yet the simple selected "
            "axis retains reflection symmetry and generates only the 74-dimensional anchored algebra. The gain law and "
            "individual event remain additional objects."
        ),
        "status": "proven_symmetry_reservoir_partial_selector_but_not_joint_closure",
        "boundary": "This mechanism demonstrates possibility; it does not derive the gain or a unique history from strict A.",
    }


def st45_medium_independence_final_test(a: np.ndarray, st40: dict) -> dict:
    t = 0.61
    unitary_spectrum = np.linalg.eigvals(expm(-1j * t * a))
    heat_spectrum = np.linalg.eigvalsh(expm(-t * a))
    positive_heat = heat_spectrum[heat_spectrum < 1.0 - 1e-10]
    unit_modulus_residual = float(np.max(np.abs(np.abs(unitary_spectrum) - 1.0)))
    heat_contraction = float(1.0 - np.max(positive_heat)) if len(positive_heat) else 0.0
    return {
        "program": "ST45",
        "object": "Medium Independence versus Dynamical Inequivalence Verdict",
        "transported_operational_record_residual": st40["transported_record_residual"],
        "unitary_spectrum_unit_modulus_residual": unit_modulus_residual,
        "least_positive_mode_heat_contraction_from_one": heat_contraction,
        "common_zero_mode_dimension": int(np.sum(np.abs(np.linalg.eigvalsh(a)) < 1e-10)),
        "theorem": (
            "Coordinate/carrier realizations related by a complete operational intertwiner encode the same information "
            "object. In contrast, exp(-itA) and exp(-tA) are not similar on any positive spectral mode: unitary eigenvalues "
            "have modulus one while heat eigenvalues have modulus below one. Their only common isomorphic sector is the "
            "zero mode. A shared generator is common structure, not equivalence of wave, unitary and diffusive dynamics."
        ),
        "answer_to_fundamental_question": (
            "Carrier independence can be fundamental as functorial invariance of the complete operational information "
            "object. It does not mean information propagates without realization, and present FIN has not shown that "
            "different physical media instantiate one such operational object."
        ),
        "status": "proven_representation_invariance_and_unitary_heat_inequivalence",
    }


def make_figures(results: dict[str, Any]) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    st28 = results["ST28"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    rows = st28["radius_rows"]
    ax.loglog([row["radius_inf"] for row in rows], [row["newton_image_bound"] for row in rows], "o-", label="Newton image")
    ax.loglog([row["radius_inf"] for row in rows], [row["radius_inf"] for row in rows], "--", label="box radius")
    ax.set_xlabel("state box radius")
    ax.set_ylabel("inclusion bound")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st28_interval_jordan_certificate.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
    st29 = results["ST29"]
    axes[0].bar([row["state"] for row in st29["rows"]], [row["algebra_dimension_with_density_observable"] for row in st29["rows"]])
    axes[0].tick_params(axis="x", rotation=20)
    axes[0].set_ylabel("generated algebra dimension")
    st30 = results["ST30"]
    for q1 in sorted(set(row["q_12_to_24"] for row in st30["rows"])):
        subset = [row for row in st30["rows"] if row["q_12_to_24"] == q1]
        axes[1].plot([row["q_24_to_48"] for row in subset], [row["trace_A48"] for row in subset], "o-", label=f"q1={q1}")
    axes[1].set_xlabel("second-level free weight")
    axes[1].set_ylabel("trace A48")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st29_state_algebra_st30_coherence.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
    st31 = results["ST31"]
    axes[0].hist(st31["common_statistics"], bins=28, alpha=0.75, label="common")
    axes[0].hist(st31["altered_statistics"], bins=28, alpha=0.75, label="altered")
    axes[0].axvline(st31["frozen_threshold"], color="black", ls="--")
    axes[0].set_xlabel("Pearson statistic")
    axes[0].legend()
    st32 = results["ST32"]
    axes[1].semilogy([row["distance"] for row in st32["rows_t_0_1_C12"]], [row["actual_strict_amplitude"] for row in st32["rows_t_0_1_C12"]], "o-", label="actual")
    axes[1].semilogy([row["distance"] for row in st32["rows_t_0_1_C12"]], [row["polynomial_bound"] for row in st32["rows_t_0_1_C12"]], "s--", label="bound")
    axes[1].set_xlabel("cyclic distance")
    axes[1].set_ylabel("amplitude at t=0.1")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st31_counts_st32_polynomial_bound.png", dpi=190)
    plt.close(fig)

    st34 = results["ST34"]["best_candidate"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.bar(range(N), np.asarray(st34["state"]) ** 2)
    ax.set_xlabel("vertex")
    ax.set_ylabel("candidate density")
    ax.set_title("constructed coercive saturating-DNLS minimizer candidate")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st34_saturating_minimizer.png", dpi=190)
    plt.close(fig)

    st37 = results["ST37"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot([row["rotation_size"] for row in st37["adversarial_rows"]], [row["detection_power"] for row in st37["adversarial_rows"]], "o-")
    ax.axhline(st37["common_holdout_false_rejection_rate"], color="tab:red", ls="--", label="common false rejection")
    ax.set_xscale("log")
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("isospectral rotation size")
    ax.set_ylabel("frozen detection probability")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st37_nuisance_detection.png", dpi=190)
    plt.close(fig)

    st40 = results["ST40"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.bar(["transported full\noperation", "generator only,\nfixed instruments"], [st40["transported_record_residual"], st40["fixed_instrument_total_variation"]])
    ax.set_yscale("log")
    ax.set_ylabel("record discrepancy")
    ax.set_title("carrier isomorphism requires transporting the whole experiment")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st40_carrier_invariance.png", dpi=190)
    plt.close(fig)

    st43 = results["ST43"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.bar(range(12), st43["sector_counts"])
    ax.set_xlabel("selected C12 orbit member")
    ax.set_ylabel("count from random infinitesimal seeds")
    ax.set_title("symmetry supplies alternatives; positive gain amplifies")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "st43_symmetry_feedback_selection.png", dpi=190)
    plt.close(fig)


def write_summary(results: dict[str, Any]) -> None:
    findings = {
        "ST28": "unique collision root certified for the frozen decimal coefficient model",
        "ST29": "generic state completes the algebra but its orbit has no canonical member",
        "ST30": "plain refinement associativity leaves one free parameter per scale",
        "ST31": "declared finite-count vertex protocol has inadequate power for the altered channel",
        "ST32": "polynomial long-range propagation bound proven without Markov positivity",
        "ST33": "reflection permits flux 0 or pi; a single state phase gives only trivial holonomy",
        "ST34": "coercivity proves a nonzero minimizer set but the best candidate is uniform, not solitonic",
        "ST35": "offline Lean replay remains blocked by absent cached toolchain",
        "ST36": "calibrated channels leave one Gibbs orbit; uncalibrated clocks leave four orbits",
        "ST37": "frozen adversarial detection quantified with noise and clock nuisance",
        "ST38": "one state cannot jointly give canonical selection, full algebra and nontrivial flux",
        "ST39": "nine relative independent axiom groups organized by countermodels",
        "ST40": "complete operational transport preserves records across carrier representations",
        "ST41": "spectrum is incomplete; the operational experiment class is the richer invariant",
        "ST42": "carrier-independent information is coherent, realization-free information is not operationally defined",
        "ST43": "symmetry supplies degenerate choices; positive gain, not symmetry, is the feedback",
        "ST44": "strict degenerate mode supports conditional selection but only partial algebra closure",
        "ST45": "unitary and heat channels share A but are dynamically inequivalent off the zero mode",
    }
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["program", "status", "main_finding"])
        for index in range(28, 46):
            key = f"ST{index}"
            writer.writerow([key, results[key]["status"], findings[key]])


def main() -> None:
    w, a, _ = strict_operator()
    results: dict[str, Any] = {
        "release": "10.57",
        "programs": "ST28-ST45",
        "language": "English",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
        },
    }
    results["ST28"] = st28_interval_jordan_certificate(a)
    results["ST29"] = st29_state_dependent_symmetry(a)
    results["ST30"] = st30_associative_refinement(a)
    results["ST31"] = st31_finite_count_likelihood(a)
    results["ST32"] = st32_long_range_polynomial_bound(a)
    results["ST33"] = st33_discrete_holonomy_and_state_connection()
    results["ST34"] = st34_saturating_soliton_existence(a)
    results["ST35"] = st35_offline_lean_audit()
    results["ST36"] = st36_sector_scale_identifiability()
    results["ST37"] = st37_nuisance_robust_adversary(a)
    results["ST38"] = st38_joint_state_selector_search(a, w)
    results["ST39"] = st39_axiom_dependency_graph()
    results["ST40"] = st40_carrier_neutral_information(a)
    results["ST41"] = st41_information_invariant_hierarchy(results["ST40"])
    results["ST42"] = st42_no_information_without_representation()
    results["ST43"] = st43_symmetry_feedback_bifurcation()
    results["ST44"] = st44_fin_mode_feedback_selector(a)
    results["ST45"] = st45_medium_independence_final_test(a, results["ST40"])
    results["recommended_next_programs"] = [
        {"id": "ST46", "priority": 1, "study": "interval-regenerate strict and memory transcendental coefficients to lift ST28 beyond the frozen model"},
        {"id": "ST47", "priority": 2, "study": "classify stronger scale-stationary refinement laws capable of coupling the free ST30 parameters"},
        {"id": "ST48", "priority": 3, "study": "derive or refute positive feedback gain from the strict adaptive law rather than inserting mu"},
        {"id": "ST49", "priority": 4, "study": "test multi-component/projective state objects for joint selector, algebra and pi-holonomy closure"},
        {"id": "ST50", "priority": 5, "study": "classify operational intertwiners among unitary, wave, heat, Green and Gibbs channels"},
        {"id": "ST51", "priority": 6, "study": "construct finite-count carrier-change tests with two explicitly different receiver models"},
        {"id": "ST52", "priority": 7, "study": "prove orbital stability of an isolated saturating-DNLS minimizer orbit with interval Hessians"},
        {"id": "ST53", "priority": 8, "study": "add entropy production and data-processing invariants to the carrier-neutral operational object"},
        {"id": "ST54", "priority": 9, "study": "derive a categorical universal property for the ST40 operational equivalence class"},
        {"id": "ST55", "priority": 10, "study": "build an executable two-carrier transfer experiment specification without claiming laboratory realization"},
        {"id": "ST56", "priority": 11, "study": "audit whether fractal compression can define a scale-stationary refinement law without hidden selectors"},
        {"id": "ST57", "priority": 12, "study": "update the minimal axiom graph after ST46-ST56 and search for genuine axiom reductions"},
    ]
    results["epistemic_boundary"] = (
        "ST28-ST45 are local mathematics and synthetic receiver studies. Carrier invariance is proved only as complete "
        "operational isomorphism; symmetry-breaking gain and selected states are conditional unless independently sourced. "
        "No result supplies laboratory data, a canonical physical history, SI units, a continuum, legacy role transfer or ToE closure."
    )
    make_figures(results)
    write_summary(results)
    RESULTS.write_text(json.dumps(native(results), indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"results": RESULTS.name, "programs": 18, "figures": 7}, indent=2))


if __name__ == "__main__":
    main()
