#!/usr/bin/env python3
"""FIN ST70--ST81: independent replay, sources, and conditional physics.

This local deterministic batch keeps exact finite theorems, rational replay,
conditional constructions, numerical continuation, synthetic receiver design,
and physical-source boundaries in separate epistemic classes.
"""

from __future__ import annotations

import base64
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
import numpy as np
import scipy
import sympy as sp
from cryptography import __version__ as cryptography_version
from cryptography.exceptions import InvalidSignature
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey, Ed25519PublicKey
from scipy.linalg import expm
from scipy.optimize import root

from fin_st01_st15_research import N, random_orthogonal_fixing_uniform, strict_operator
from fin_st28_st45_research import dyadic_lift, saturation_energy_gradient_hessians
from fin_st46_st57_research import carrier_probability_table


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST70_ST81_Results.json"
SUMMARY = ROOT / "FIN_ST70_ST81_Summary.csv"
REPLAY70 = ROOT / "FIN_ST70_Rational_Replay_Certificate.json"
BOUNDS74 = ROOT / "FIN_ST74_Robust_Count_Design.json"
SPEC80 = ROOT / "FIN_ST80_Signed_Custody_Validator.json"
FIG_DIR = ROOT / "FIN_ST70_ST81_Figures"
SEED = 20260815


def native(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(item) for item in value]
    if isinstance(value, np.ndarray):
        return native(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, Fraction):
        return str(value)
    return value


def canonical_bytes(payload: Any) -> bytes:
    return json.dumps(native(payload), sort_keys=True, separators=(",", ":")).encode("utf-8")


def canonical_digest(payload: Any) -> str:
    return hashlib.sha256(canonical_bytes(payload)).hexdigest()


def frac(value: Any) -> Fraction:
    return Fraction(str(value))


def fi_add(*values: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return sum((value[0] for value in values), Fraction(0)), sum(
        (value[1] for value in values), Fraction(0)
    )


def fi_mul(
    left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    products = [left[i] * right[j] for i in range(2) for j in range(2)]
    return min(products), max(products)


def st70_exact_rational_replay() -> dict:
    """Replay every decisive exported ST58 inequality with exact Fractions.

    This deliberately does not call mpmath, numpy interval helpers, or the
    original ST58 implementation.  It verifies the certificate ledger and the
    final root bracket exactly as decimal rationals, while retaining the
    upstream transcendental-enclosure trust boundary.
    """
    certificate = json.loads((ROOT / "FIN_ST58_Full_Interval_Certificate.json").read_text(encoding="utf-8"))
    coefficients = certificate["coefficient_intervals"]
    a = tuple(map(frac, coefficients["a=c2+c3"]))
    b = tuple(map(frac, coefficients["b=2c1"]))
    c = tuple(map(frac, coefficients["c=c0"]))
    speed = tuple(map(frac, certificate["collision_speed_interval"]))

    def polynomial_at(point: Fraction) -> tuple[Fraction, Fraction]:
        scalar = (point, point)
        return fi_add(a, fi_mul(b, scalar), fi_mul(c, fi_mul(scalar, scalar)))

    left_sign = polynomial_at(speed[0])
    right_sign = polynomial_at(speed[1])
    derivative = fi_add(b, fi_mul((2 * c[0], 2 * c[1]), speed))
    exact_checks = {
        "left_polynomial_strictly_positive": left_sign[0] > 0,
        "right_polynomial_strictly_negative": right_sign[1] < 0,
        "derivative_strictly_negative_on_interval": derivative[1] < 0,
        "stationary_krawczyk_margin_positive": frac(certificate["stationary_root"]["minimum_inclusion_margin"]) > 0,
        "stationary_defect_below_one": frac(certificate["stationary_root"]["defect_infinity_norm_upper"]) < 1,
        "linear_inclusion_margins_positive": all(frac(value) > 0 for value in certificate["linear_solve_minimum_inclusion_margins"]),
        "linear_comparison_radii_below_one": all(frac(value) < 1 for value in certificate["linear_solve_comparison_spectral_radii"]),
        "poles_positive": all(frac(interval[0]) > 0 for interval in certificate["pole_intervals"]),
        "residues_positive": all(frac(interval[0]) > 0 for interval in certificate["residue_intervals"]),
    }
    packet = {
        "input_sha256": hashlib.sha256((ROOT / "FIN_ST58_Full_Interval_Certificate.json").read_bytes()).hexdigest(),
        "arithmetic": "Python standard-library Fraction; every serialized decimal is interpreted exactly",
        "root_equation": "p(s)=a+b*s+c*s^2 after x=1/s",
        "left_polynomial_interval": [str(left_sign[0]), str(left_sign[1])],
        "right_polynomial_interval": [str(right_sign[0]), str(right_sign[1])],
        "derivative_interval": [str(derivative[0]), str(derivative[1])],
        "checks": exact_checks,
        "all_checks_pass": all(exact_checks.values()),
        "collision_speed_interval": certificate["collision_speed_interval"],
        "trust_boundary": (
            "The replay is independent of the ST58 interval code and exactly verifies its exported acceptance inequalities. "
            "It does not independently regenerate cosine, fifth-root, Fourier, matrix-preconditioner, or stationary-box enclosures."
        ),
    }
    REPLAY70.write_text(json.dumps(packet, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST70",
        "object": "Standard-Library Exact-Rational Replay of the ST58 Certificate",
        **packet,
        "status": "proven_exact_rational_replay_with_explicit_upstream_interval_trust_boundary",
    }


def reflection_odd_basis() -> np.ndarray:
    basis = []
    for site in range(1, 6):
        vector = np.zeros(N)
        vector[site] = 1.0 / math.sqrt(2.0)
        vector[-site] = -1.0 / math.sqrt(2.0)
        basis.append(vector)
    return np.column_stack(basis)


def st71_time_oriented_response_classification(a: np.ndarray) -> dict:
    odd = reflection_odd_basis()
    reflection = np.zeros((N, N))
    for x in range(N):
        reflection[(-x) % N, x] = 1.0
    odd_residual = float(np.linalg.norm(reflection @ odd + odd))
    reduced_a = odd.T @ a @ odd
    rows = []
    for curvature in [-2.0, -0.25, 0.0, 0.25, 2.0]:
        response_hessian = curvature * np.eye(odd.shape[1])
        gain_generator = -response_hessian
        rows.append(
            {
                "odd_sector_curvature": curvature,
                "gradient_flow_gain": -curvature,
                "uniform_positive_gain": bool(np.min(np.linalg.eigvalsh(gain_generator)) > 0),
                "oriented_response_inequality_with_mu_0_2": bool(curvature <= -0.2),
            }
        )
    return {
        "program": "ST71",
        "object": "Time-Oriented Odd-Response Cone Classification",
        "reflection_odd_state_dimension": odd.shape[1],
        "reflection_odd_basis_residual": odd_residual,
        "strict_odd_sector_eigenvalues": np.linalg.eigvalsh(reduced_a),
        "rows": rows,
        "theorem": (
            "Let xi be the five-dimensional reflection-odd perturbation and let H_- be the Hessian of an adaptive "
            "response at the stationary strict operator. The linearized gradient law is dot(xi)=-H_- xi. Uniform "
            "positive gain at rate at least mu is therefore equivalent to H_- <= -mu I. Removing this oriented "
            "one-sided inequality admits H_- and -H_- and restores the ST59 sign ambiguity."
        ),
        "necessity": "A one-sided time/order premise is necessary and sufficient for uniform gain within this declared linear response class.",
        "status": "proven_relative_classification_of_time_oriented_positive_gain_cone",
        "boundary": "The inequality H_-<=-mu I is an added response law; strict A and its reflection symmetry do not source its sign or mu.",
    }


def gf2_rank(matrix: np.ndarray) -> int:
    work = np.asarray(matrix, dtype=np.uint8).copy() % 2
    row = 0
    for col in range(work.shape[1]):
        pivots = np.flatnonzero(work[row:, col])
        if len(pivots) == 0:
            continue
        pivot = row + int(pivots[0])
        work[[row, pivot]] = work[[pivot, row]]
        for other in range(work.shape[0]):
            if other != row and work[other, col]:
                work[other] ^= work[row]
        row += 1
        if row == work.shape[0]:
            break
    return row


def st72_spin_central_extension_obstruction(a: np.ndarray) -> dict:
    n = 12

    def cocycle(left: int, right: int) -> int:
        return int(left + right >= n)

    cocycle_defects = 0
    for left in range(n):
        for middle in range(n):
            for right in range(n):
                defect = (
                    cocycle(left, middle)
                    + cocycle((left + middle) % n, right)
                    - cocycle(middle, right)
                    - cocycle(left, (middle + right) % n)
                ) % 2
                cocycle_defects += int(defect != 0)

    # Test c=delta f over every normalized one-cochain f(0)=0 by GF(2) ranks.
    equations, target = [], []
    for left in range(n):
        for right in range(n):
            row = np.zeros(n - 1, dtype=np.uint8)
            for index in [left, right, (left + right) % n]:
                if index != 0:
                    row[index - 1] ^= 1
            equations.append(row)
            target.append(cocycle(left, right))
    equations = np.asarray(equations, dtype=np.uint8)
    target = np.asarray(target, dtype=np.uint8)[:, None]
    rank = gf2_rank(equations)
    augmented_rank = gf2_rank(np.hstack([equations, target]))

    # A scalar representation and its lift have the same action on the base.
    base_rotation = np.roll(np.eye(N), 1, axis=0)
    scalar_commutator = float(np.linalg.norm(a @ base_rotation - base_rotation @ a))
    return {
        "program": "ST72",
        "object": "C24 Spin-Cover Construction and Strict-Source Obstruction",
        "normalized_cocycle": "c(a,b)=floor((a+b)/12) mod 2",
        "cocycle_equations_tested": n**3,
        "cocycle_defects": cocycle_defects,
        "coboundary_matrix_rank": rank,
        "coboundary_augmented_rank": augmented_rank,
        "cocycle_is_nontrivial": bool(augmented_rank > rank),
        "extension_classes": ["C12 x Z2 (lift^12=+1)", "C24 (lift^12=-1)"],
        "strict_rotation_commutator_residual": scalar_commutator,
        "theorem": (
            "A central Z2 extension of cyclic C12 is classified by the sign of the twelfth power of a lifted generator. "
            "Because 12 is even, changing the lift by the central sign cannot change that invariant. The displayed carry "
            "cocycle is exact-cocycle verified and is not a coboundary, so the nontrivial C24 spin cover exists. Both "
            "extension classes project to the same scalar C12 action; a scalar function of strict A cannot select one."
        ),
        "status": "proven_nontrivial_spin_cover_exists_and_scalar_strict_data_do_not_select_it",
        "boundary": "The C24 lift is a mathematically admissible added object, not a strict-derived spin structure or QW-2191 discharge.",
    }


def st73_nonstationary_fine_data_refinement(a: np.ndarray) -> dict:
    n = len(a)
    fine = np.vstack([np.eye(n), -np.eye(n)]) / math.sqrt(2.0)
    baseline = float(np.trace(fine.T @ dyadic_lift(a, 0.0) @ fine) / n)
    slope = float(
        (np.trace(fine.T @ dyadic_lift(a, 1.0) @ fine) / n) - baseline
    )
    penalty = 0.6
    rows = []
    for true_q in [0.0, 0.15, 0.4, 0.9]:
        fine_statistic = slope * true_q
        selected = max(0.0, slope * fine_statistic / (slope**2 + penalty))
        rows.append(
            {
                "supplied_fine_q": true_q,
                "fine_only_trace_increment": fine_statistic,
                "rate_distortion_selected_q": selected,
                "shrinkage_bias": selected - true_q,
            }
        )

    # Online nonstationary gradient descent for a frozen illustrative stream.
    stream = [0.10] * 20 + [0.75] * 20 + [0.25] * 20
    q = 0.0
    learning_rate = 0.08
    trajectory = []
    for time, supplied_q in enumerate(stream):
        statistic = slope * supplied_q
        gradient = slope * (slope * q - statistic) + penalty * q
        q = max(0.0, q - learning_rate * gradient)
        trajectory.append({"time": time, "supplied_q": supplied_q, "selected_q": q})
    return {
        "program": "ST73",
        "object": "Nonstationary Fine-Mode Rate-Distortion Refinement",
        "fine_subspace_dimension": n,
        "fine_trace_baseline": baseline,
        "fine_statistic_slope_in_q": slope,
        "quadratic_rate_penalty": penalty,
        "closed_form_selector": "q* = max(0, g*y/(g^2+lambda))",
        "rows": rows,
        "online_trajectory": trajectory,
        "theorem": (
            "Once a fine-only statistic y and penalty lambda>0 are supplied, D(q)=|gq-y|^2/2+lambda q^2/2 is "
            "strictly convex and uniquely selects q*=max(0,g y/(g^2+lambda)). Conversely every q0>=0 is selected "
            "by y=(g^2+lambda)q0/g. Thus nonstationary fine data can remove the lift degeneracy, but the answer is "
            "carried by the added fine record and coding penalty rather than coarse strict A."
        ),
        "status": "proven_conditional_fine_data_selection_and_no_canonical_source_from_coarse_operator",
        "boundary": "The fine statistic, penalty and data stream are constructed; no FIN-internal fractal code or canonical q is derived.",
    }


def st74_nuisance_robust_count_design(a: np.ndarray) -> dict:
    rng = np.random.default_rng(20260813 + 51)
    transport = random_orthogonal_fixing_uniform(rng)
    p = carrier_probability_table(a, np.eye(N), transported=False) / N
    q = carrier_probability_table(a, transport, transported=False) / N
    uniform = np.full_like(p, 1.0 / (N * N))
    maximum_loss_plus_dark = 0.12
    grid = np.linspace(0.0, maximum_loss_plus_dark, 241)
    best: dict[str, float] | None = None
    for epsilon_p in grid:
        mixed_p = (1.0 - epsilon_p) * p + epsilon_p * uniform
        for epsilon_q in grid:
            mixed_q = (1.0 - epsilon_q) * q + epsilon_q * uniform
            l1 = float(np.sum(np.abs(mixed_p - mixed_q)))
            if best is None or l1 < best["l1"]:
                kl_pq = float(np.sum(mixed_p * np.log(mixed_p / mixed_q)))
                kl_qp = float(np.sum(mixed_q * np.log(mixed_q / mixed_p)))
                best = {
                    "l1": l1,
                    "epsilon_p": float(epsilon_p),
                    "epsilon_q": float(epsilon_q),
                    "kl_pq": kl_pq,
                    "kl_qp": kl_qp,
                }
    assert best is not None
    calibration_tv_each = 0.003
    robust_l1 = max(0.0, best["l1"] - 4.0 * calibration_tv_each)
    target_error = 0.01
    alphabet = N * N
    total_sufficient = math.ceil(
        8.0 * (alphabet * math.log(2.0) + math.log(1.0 / target_error)) / robust_l1**2
    )
    binary_kl = (1.0 - 2.0 * target_error) * math.log((1.0 - target_error) / target_error)
    total_necessary = math.ceil(max(binary_kl / best["kl_pq"], binary_kl / best["kl_qp"]))
    packet = {
        "synthetic_table_source": "ST51 deterministic carrier tables",
        "loss_plus_dark_uniform_mixture_range": [0.0, maximum_loss_plus_dark],
        "nuisance_grid_points_per_hypothesis": len(grid),
        "calibration_TV_ball_each_hypothesis": calibration_tv_each,
        "closest_grid_pair": best,
        "robust_joint_L1_separation": robust_l1,
        "joint_alphabet_size": alphabet,
        "target_each_error": target_error,
        "necessary_total_shots_for_selected_adversarial_pair": total_necessary,
        "necessary_mean_shots_per_preparation_ceiling": math.ceil(total_necessary / N),
        "distribution_free_sufficient_total_shots": total_sufficient,
        "distribution_free_sufficient_mean_shots_per_preparation_ceiling": math.ceil(total_sufficient / N),
        "sampling_design": "each joint shot chooses one of 12 preparations uniformly; per-preparation values are expected means, not fixed-stratum guarantees",
        "method": "KL data-processing lower bound plus multinomial L1 concentration upper bound",
    }
    BOUNDS74.write_text(json.dumps(packet, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST74",
        "object": "Nuisance-Robust Synthetic Minimax Count Bracket",
        **packet,
        "status": "proven_conservative_count_bracket_for_declared_synthetic_uncertainty_sets",
        "boundary": "The nuisance geometry is an assumed uniform-mixture/TV model, not detector calibration or a laboratory count recommendation.",
    }


def spectral_groups(a: np.ndarray) -> tuple[np.ndarray, list[np.ndarray], list[int]]:
    eigenvalues, vectors = np.linalg.eigh(a)
    groups: list[list[int]] = []
    for index, value in enumerate(eigenvalues):
        if not groups or abs(value - eigenvalues[groups[-1][0]]) > 1e-9:
            groups.append([index])
        else:
            groups[-1].append(index)
    projectors = [vectors[:, group] @ vectors[:, group].T for group in groups]
    return eigenvalues, projectors, [len(group) for group in groups]


def st75_all_cptp_intertwiners(a: np.ndarray) -> dict:
    eigenvalues, projectors, multiplicities = spectral_groups(a)
    block_dimension = sum(size**2 for size in multiplicities)
    t = tau = 0.61
    maximum_phase = t * float(eigenvalues[-1] - eigenvalues[0])
    pinching = lambda rho: sum(projector @ rho @ projector for projector in projectors)
    sigma = expm(-0.8 * a)
    sigma /= np.trace(sigma)
    reset = lambda rho: np.trace(rho) * sigma

    rng = np.random.default_rng(SEED + 75)
    trials = []
    unitary = expm(-1j * t * a)

    def dephase(rho: np.ndarray) -> np.ndarray:
        values, vectors = np.linalg.eigh(a)
        spectral = vectors.T.conj() @ rho @ vectors
        gaps = values[:, None] - values[None, :]
        return vectors @ (np.exp(-tau * gaps**2) * spectral) @ vectors.T.conj()

    for name, channel in [("pinching", pinching), ("Gibbs reset", reset)]:
        residuals = []
        for _ in range(8):
            matrix = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
            rho = matrix @ matrix.conj().T
            rho /= np.trace(rho)
            left = channel(unitary @ rho @ unitary.conj().T)
            right = dephase(channel(rho))
            residuals.append(float(np.linalg.norm(left - right)))
        trials.append({"channel": name, "maximum_intertwining_residual": max(residuals)})
    return {
        "program": "ST75",
        "object": "Complete CPTP Intertwiner Cone Characterization",
        "parameters": {"unitary_time": t, "dephasing_time": tau},
        "maximum_nonzero_unitary_gap_phase": maximum_phase,
        "below_two_pi": bool(maximum_phase < 2.0 * math.pi),
        "spectral_block_multiplicities": multiplicities,
        "fixed_block_algebra_dimension": block_dimension,
        "complex_linear_intertwiner_space_dimension": block_dimension**2,
        "choi_semidefinite_characterization": "J(Phi)>=0, Tr_out J(Phi)=I, Phi=Pi o Phi o Pi",
        "example_channels": trials,
        "theorem": (
            "For Phi o U_t=D_tau o Phi, every input eigenoperator has a unit-modulus phase while every nonfixed "
            "dephasing eigenoperator has modulus below one. Since the strict phase range is below 2pi, equality is "
            "possible only for zero input and output gaps. Hence every linear intertwiner is exactly a map from the "
            "22-dimensional equal-energy block algebra to itself, extended by zero on cross-energy coherences. The "
            "CPTP intertwiners are precisely the Choi-positive trace-preserving members of this 484-dimensional linear space."
        ),
        "status": "proven_complete_finite_CPTP_intertwiner_cone_characterization",
        "boundary": "This classifies unitary versus A-dephasing, not the classical vertex heat semigroup and not physical carrier equivalence.",
    }


def st76_thermodynamic_resource_minimality() -> dict:
    rows = [
        {
            "resource": "positive energy scale E_*",
            "removal_countermodel": "H and cH give identical dimensionless A records while assigning different work/energy",
            "necessary_for": "dimensionful Hamiltonian and work",
        },
        {
            "resource": "entropy unit k_*",
            "removal_countermodel": "S_vN and c S_vN are indistinguishable as dimensionless information",
            "necessary_for": "physical entropy and temperature conversion",
        },
        {
            "resource": "full-rank Gibbs reference/bath state gamma",
            "removal_countermodel": "without gamma no preferred beta or free-energy zero is selected",
            "necessary_for": "equilibrium comparison and heat convention",
        },
        {
            "resource": "gamma-preserving process class and heat/work instrument",
            "removal_countermodel": "an arbitrary channel may prepare a higher-resource state from gamma",
            "necessary_for": "monotonicity and operational entropy production",
        },
    ]
    return {
        "program": "ST76",
        "object": "Relative Minimality of the Thermodynamic Conversion Package",
        "resources": rows,
        "resource_count": len(rows),
        "scale_gauge": "(E_*,T)->(c E_*,c T) with beta_tilde fixed",
        "theorem": (
            "Relative to the target of dimensionful energy, entropy, Gibbs equilibrium and nonnegative operational "
            "entropy production, the four displayed resource types are independently necessary. Their joint presence "
            "is sufficient for the standard identity D(rho||gamma)=beta[F(rho)-F(gamma)] and its data-processing "
            "monotonicity under gamma-preserving channels. Removing each item admits the stated countermodel."
        ),
        "status": "proven_relative_resource_minimality_with_removal_countermodels",
        "boundary": "The theorem organizes added conversion resources; it does not derive E_*, k_B, a bath, or a physical instrument from FIN.",
    }


def reflection_expansion() -> np.ndarray:
    expansion = np.zeros((N, 7))
    expansion[0, 0] = 1.0
    for site in range(1, 6):
        expansion[site, site] = 1.0
        expansion[-site, site] = 1.0
    expansion[6, 6] = 1.0
    return expansion


def st77_pseudo_arclength_fold(a: np.ndarray) -> dict:
    expansion = reflection_expansion()
    selected_rows = np.arange(7)

    def reduced(q: np.ndarray, kappa: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        state = expansion @ q
        _, gradient, hessian, _ = saturation_energy_gradient_hessians(kappa * a, state)
        jacobian = hessian[np.ix_(selected_rows, np.arange(N))] @ expansion
        parameter_column = (a @ state)[selected_rows]
        return gradient[selected_rows], jacobian, parameter_column

    q = np.zeros(7)
    q[0] = 2.669532756662446
    q = root(lambda value: reduced(value, 0.0)[0], q, jac=lambda value: reduced(value, 0.0)[1], options={"xtol": 1e-12}).x
    _, jacobian, parameter_column = reduced(q, 0.0)
    tangent = np.r_[np.linalg.solve(jacobian, -parameter_column), 1.0]
    tangent /= np.linalg.norm(tangent)
    point = np.r_[q, 0.0]
    step_size = 5e-4
    rows = []
    passed_fold = False
    fold_seed: np.ndarray | None = None
    fold_tangent: np.ndarray | None = None
    previous_tangent = tangent.copy()

    for step in range(2000):
        predictor = point + step_size * tangent

        def augmented(candidate: np.ndarray) -> np.ndarray:
            gradient = reduced(candidate[:-1], candidate[-1])[0]
            return np.r_[gradient, np.dot(candidate - predictor, tangent)]

        def augmented_jacobian(candidate: np.ndarray) -> np.ndarray:
            _, state_jacobian, column = reduced(candidate[:-1], candidate[-1])
            return np.block([[state_jacobian, column[:, None]], [tangent[None, :]]])

        solution = root(augmented, predictor, jac=augmented_jacobian, method="lm", options={"ftol": 1e-13, "xtol": 1e-13, "gtol": 1e-13, "maxiter": 1000})
        residual = float(np.linalg.norm(augmented(solution.x), ord=np.inf))
        if residual > 1e-9:
            step_size *= 0.5
            if step_size < 1e-7:
                break
            continue
        new_point = solution.x
        _, state_jacobian, column = reduced(new_point[:-1], new_point[-1])
        _, singular_values, vh = np.linalg.svd(np.column_stack([state_jacobian, column]))
        new_tangent = vh[-1]
        if np.dot(new_tangent, tangent) < 0:
            new_tangent = -new_tangent
        new_tangent /= np.linalg.norm(new_tangent)
        state = expansion @ new_point[:-1]
        power = float(state @ state)
        ipr = float(np.sum(state**4) / power**2)
        if step % 4 == 0:
            rows.append(
                {
                    "arclength_step": step,
                    "kappa": float(new_point[-1]),
                    "IPR": ipr,
                    "tangent_kappa_component": float(new_tangent[-1]),
                    "reduced_J_minimum_singular_value": float(np.linalg.svd(state_jacobian, compute_uv=False)[-1]),
                    "residual_inf": residual,
                }
            )
        if previous_tangent[-1] > 0.0 and new_tangent[-1] <= 0.0 and not passed_fold:
            passed_fold = True
            fold_seed = new_point.copy()
            fold_tangent = new_tangent.copy()
        point, previous_tangent, tangent = new_point, new_tangent, new_tangent
        if passed_fold and point[-1] < 0.012:
            break

    if fold_seed is None or fold_tangent is None:
        raise RuntimeError("ST77 did not traverse the first fold")
    q_seed, kappa_seed = fold_seed[:-1], float(fold_seed[-1])
    _, fold_jacobian, _ = reduced(q_seed, kappa_seed)
    null_vector = np.linalg.svd(fold_jacobian)[2][-1]
    null_vector /= np.linalg.norm(null_vector)
    fold_start = np.r_[q_seed, kappa_seed, null_vector]

    def fold_system(candidate: np.ndarray) -> np.ndarray:
        q_value, kappa_value, vector = candidate[:7], float(candidate[7]), candidate[8:]
        gradient, matrix, _ = reduced(q_value, kappa_value)
        return np.r_[gradient, matrix @ vector, 0.5 * (vector @ vector - 1.0)]

    fold_solution = root(fold_system, fold_start, method="lm", options={"ftol": 1e-14, "xtol": 1e-14, "gtol": 1e-14, "maxiter": 10000})
    fold_residual = float(np.linalg.norm(fold_system(fold_solution.x), ord=np.inf))
    fold_state = expansion @ fold_solution.x[:7]
    fold_power = float(fold_state @ fold_state)

    # Numerical transversality of the augmented fold equations.
    dimension = len(fold_solution.x)
    numerical_jacobian = np.empty((dimension, dimension))
    for column_index in range(dimension):
        delta = 1e-6 * max(1.0, abs(float(fold_solution.x[column_index])))
        direction = np.zeros(dimension)
        direction[column_index] = delta
        numerical_jacobian[:, column_index] = (
            fold_system(fold_solution.x + direction) - fold_system(fold_solution.x - direction)
        ) / (2.0 * delta)
    augmented_min_sv = float(np.linalg.svd(numerical_jacobian, compute_uv=False)[-1])
    return {
        "program": "ST77",
        "object": "Reflection-Reduced Pseudo-Arclength Continuation Through the ST65 Fold",
        "stored_branch_rows": rows,
        "fold_kappa": float(fold_solution.x[7]),
        "fold_state_reduced": fold_solution.x[:7],
        "fold_null_vector": fold_solution.x[8:],
        "fold_IPR": float(np.sum(fold_state**4) / fold_power**2),
        "fold_system_residual_inf": fold_residual,
        "augmented_fold_jacobian_minimum_singular_value": augmented_min_sv,
        "return_branch_last_kappa": float(point[-1]),
        "return_branch_found": bool(passed_fold and point[-1] < 0.012),
        "status": "strong_numerical_simple_fold_and_return_branch_not_interval_certified",
        "boundary": "Pseudo-arclength corrects the fixed-kappa stopping interpretation but does not interval-certify the fold, reach kappa=1, or exclude disconnected branches.",
    }


def st78_strict_doublet_backreaction(a: np.ndarray) -> dict:
    eigenvalues = np.linalg.eigvalsh(a)
    lambda_one = float(eigenvalues[1])
    mu = 7.0 / 20.0
    epsilon = 1.0 / 100.0
    nu = 1.0 / 100.0
    critical_coupling = mu / lambda_one

    def radial_polynomial(y: float, coupling: float) -> float:
        return -(mu - coupling * lambda_one) + y - epsilon * y**5 + nu * y**6

    rows = []
    for coupling in [0.0, 0.25, 0.999 * critical_coupling, 1.001 * critical_coupling, 1.0]:
        effective_mu = mu - coupling * lambda_one
        if effective_mu > 0:
            solution = root(lambda value: [radial_polynomial(float(value[0]), coupling)], [effective_mu]).x[0]
            root_count = 1
            positive_root = float(solution)
        else:
            root_count = 0
            positive_root = None
        rows.append(
            {
                "coupling": coupling,
                "effective_mu": effective_mu,
                "positive_radial_root_count": root_count,
                "positive_radial_root": positive_root,
            }
        )
    return {
        "program": "ST78",
        "object": "Strict-Doublet Stiffness Backreaction on the Twelve-Branch Potential",
        "strict_first_doublet_eigenvalue": lambda_one,
        "critical_coupling_mu_over_lambda1": critical_coupling,
        "rows": rows,
        "theorem": (
            "Embedding the ST66 order parameter in the first strict Fourier doublet adds g lambda_1 |z|^2/2 and "
            "replaces mu by mu-g lambda_1. The remaining radial function y-0.01 y^5+0.01 y^6 is strictly increasing "
            "and nonnegative for y>=0. Therefore exactly one nonzero radial branch exists for g<mu/lambda_1, and none "
            "exists for g>=mu/lambda_1. At unit coupling strict stiffness quenches the constructed twelve-branch phase."
        ),
        "status": "proven_conditional_backreaction_threshold_for_constructed_potential",
        "boundary": "The coupling and potential are added structures; quenching is not a theorem about every FIN nonlinearity or physical phase.",
    }


def bargmann_sum(rays: np.ndarray) -> float:
    value = 0.0
    for x in range(len(rays)):
        left, middle, right = rays[x], rays[(x + 1) % len(rays)], rays[(x + 2) % len(rays)]
        value += float(np.imag(np.vdot(left, middle) * np.vdot(middle, right) * np.vdot(right, left)))
    return value


def st79_orientation_odd_bargmann_obstruction(a: np.ndarray) -> dict:
    reflection = (-np.arange(N)) % N
    rows = []
    for time in [0.1, 0.63, 1.0, 2.0]:
        rays = expm(-1j * time * a)
        chirality = bargmann_sum(rays)
        reflected = bargmann_sum(rays[reflection])
        rows.append(
            {
                "time": time,
                "strict_functional_calculus_chirality": chirality,
                "reflected_chirality": reflected,
                "odd_sum_residual": abs(chirality + reflected),
            }
        )
    return {
        "program": "ST79",
        "object": "Stabilizer No-Go for a Strict-Sourced Bargmann Chirality",
        "rows": rows,
        "maximum_strict_candidate_absolute_chirality": max(abs(row["strict_functional_calculus_chirality"]) for row in rows),
        "theorem": (
            "Reflection belongs to the stabilizer of real circulant strict A. Any deterministic stabilizer-covariant "
            "ray construction has the same output orbit under reflection, while the imaginary Bargmann sum changes "
            "sign. Its value must therefore equal its negative and vanish. The sampled functional-calculus candidates "
            "confirm the exact symmetry conclusion to floating residual."
        ),
        "status": "proven_stabilizer_obstruction_for_deterministic_strict_bargmann_chirality",
        "boundary": "A state, boundary, projective lift or symmetry-breaking record can evade the theorem only as an additional sourced object.",
    }


def raw_public_key(private_key: Ed25519PrivateKey) -> bytes:
    return private_key.public_key().public_bytes(
        encoding=serialization.Encoding.Raw, format=serialization.PublicFormat.Raw
    )


def st80_signed_custody_validator(st74: dict) -> dict:
    roles = ["provider", "registrar", "analyst", "custodian"]
    private_keys = {
        role: Ed25519PrivateKey.from_private_bytes(hashlib.sha256(f"FIN-ST80-test-key:{role}".encode()).digest())
        for role in roles
    }
    public_keys = {role: raw_public_key(private_keys[role]) for role in roles}
    public_b64 = {role: base64.b64encode(public_keys[role]).decode("ascii") for role in roles}
    fingerprints = {role: hashlib.sha256(public_keys[role]).hexdigest() for role in roles}
    events = [
        {
            "timestamp": f"2026-08-11T12:00:{index:02d}Z",
            "carrier_id": "synthetic-A",
            "preparation_id": index % N,
            "outcome_id": (5 * index + 1) % N,
            "run_id": "ST80-test-run",
            "calibration_split": "holdout" if index >= 8 else "calibration",
            "blind_label": f"B{index:02d}",
        }
        for index in range(12)
    ]
    nuisance_payload = {
        "loss_plus_dark_range": st74["loss_plus_dark_uniform_mixture_range"],
        "calibration_TV_ball": st74["calibration_TV_ball_each_hypothesis"],
        "target_error": st74["target_each_error"],
    }
    core = {
        "version": "1.0.0",
        "model": "FIN ST74 synthetic nuisance-count design",
        "events_sha256": canonical_digest(events),
        "nuisance_sha256": canonical_digest(nuisance_payload),
        "public_keys": public_b64,
        "role_order": roles,
        "holdout_frozen": True,
        "single_pipeline_run": True,
        "no_refit_after_failure": True,
    }
    core_hash = canonical_digest(core)
    attestations = []
    previous_hash = "0" * 64
    for stage, role in enumerate(roles):
        message = {
            "stage": stage,
            "role": role,
            "core_hash": core_hash,
            "previous_attestation_hash": previous_hash,
            "public_key_fingerprint": fingerprints[role],
        }
        signature = private_keys[role].sign(canonical_bytes(message))
        attestation = {
            "message": message,
            "signature_base64": base64.b64encode(signature).decode("ascii"),
        }
        attestations.append(attestation)
        previous_hash = canonical_digest(attestation)

    packet = {"core": core, "events": events, "nuisance": nuisance_payload, "attestations": attestations}

    def validate(candidate: dict) -> dict:
        checks: dict[str, bool] = {}
        candidate_core = candidate["core"]
        candidate_events = candidate["events"]
        candidate_nuisance = candidate["nuisance"]
        checks["event_hash"] = canonical_digest(candidate_events) == candidate_core["events_sha256"]
        checks["nuisance_hash"] = canonical_digest(candidate_nuisance) == candidate_core["nuisance_sha256"]
        key_bytes = {
            role: base64.b64decode(candidate_core["public_keys"][role])
            for role in candidate_core["role_order"]
        }
        checks["distinct_public_keys"] = len(set(key_bytes.values())) == len(roles)
        expected_previous = "0" * 64
        signature_checks = []
        chain_checks = []
        for stage, attestation in enumerate(candidate["attestations"]):
            message = attestation["message"]
            role = message["role"]
            chain_checks.append(
                message["stage"] == stage
                and role == roles[stage]
                and message["previous_attestation_hash"] == expected_previous
                and message["core_hash"] == canonical_digest(candidate_core)
                and message["public_key_fingerprint"] == hashlib.sha256(key_bytes[role]).hexdigest()
            )
            try:
                Ed25519PublicKey.from_public_bytes(key_bytes[role]).verify(
                    base64.b64decode(attestation["signature_base64"]), canonical_bytes(message)
                )
                signature_checks.append(True)
            except (InvalidSignature, ValueError):
                signature_checks.append(False)
            expected_previous = canonical_digest(attestation)
        checks["attestation_chain"] = all(chain_checks) and len(chain_checks) == len(roles)
        checks["signatures"] = all(signature_checks) and len(signature_checks) == len(roles)
        checks["raw_events_present"] = len(candidate_events) > 0
        checks["frozen_protocol"] = bool(candidate_core["holdout_frozen"] and candidate_core["single_pipeline_run"] and candidate_core["no_refit_after_failure"])
        return {"checks": checks, "accepted": all(checks.values())}

    valid = validate(json.loads(json.dumps(packet)))
    tampered = json.loads(json.dumps(packet))
    tampered["events"][0]["outcome_id"] = (tampered["events"][0]["outcome_id"] + 1) % N
    tampered_result = validate(tampered)
    duplicate = json.loads(json.dumps(packet))
    duplicate["core"]["public_keys"]["registrar"] = duplicate["core"]["public_keys"]["provider"]
    duplicate_result = validate(duplicate)
    output = {
        "specification": packet,
        "packet_sha256": canonical_digest(packet),
        "valid_case": valid,
        "tampered_event_case": tampered_result,
        "duplicate_key_case": duplicate_result,
        "test_key_warning": "Keys are deterministic local test fixtures and provide no real-world identity or independent custody.",
    }
    SPEC80.write_text(json.dumps(output, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST80",
        "object": "Ed25519 Signed-Custody and Frozen-Nuisance Validator",
        "packet_sha256": output["packet_sha256"],
        "role_public_key_fingerprints": fingerprints,
        "valid_case_accepted": valid["accepted"],
        "tampered_event_rejected": not tampered_result["accepted"],
        "duplicate_public_key_rejected": not duplicate_result["accepted"],
        "attestation_count": len(attestations),
        "status": "constructed_executable_signed_custody_validator_with_synthetic_test_keys",
        "boundary": output["test_key_warning"] + " The validator cannot create an external laboratory event or trusted human identity.",
    }


def st81_architecture_reconciliation() -> dict:
    rows = [
        {"source_group": "A0 strict finite operator", "package": "W0", "post_batch": "retained", "evidence": "ST70 audits consequences, not provenance"},
        {"source_group": "A1 state-selection event", "package": "SA", "post_batch": "retained", "evidence": "ST72/ST79 obstruct strict spin/chirality selection"},
        {"source_group": "A2 refinement law/fine record", "package": "RA", "post_batch": "retained", "evidence": "ST73 selects only after supplied fine data and penalty"},
        {"source_group": "A3 dimensional calibration", "package": "CA", "post_batch": "retained", "evidence": "ST76 proves relative conversion-resource necessity"},
        {"source_group": "A4 operational process", "package": "OA", "post_batch": "constructible", "evidence": "ST75 characterizes one CPTP cone; ST80 validates syntax"},
        {"source_group": "A5 external custody/data", "package": "OA-external", "post_batch": "retained", "evidence": "ST80 test keys are local and synthetic"},
        {"source_group": "A6 nonlinear response/gain", "package": "SA", "post_batch": "retained", "evidence": "ST71 needs oriented curvature; ST78 coupling is inserted"},
        {"source_group": "A7 connection/projective/spin sector", "package": "SA", "post_batch": "retained", "evidence": "ST72 constructs two covers but strict scalar data do not choose"},
        {"source_group": "A8 temporal memory realization", "package": "SA/TA", "post_batch": "retained", "evidence": "ST71 adds time order but does not derive a temporal realization"},
    ]
    return {
        "program": "ST81",
        "object": "Post-ST80 W0+CA+SA Source Architecture",
        "rows": rows,
        "strict_source_group_count_before": 9,
        "strict_source_group_count_after": 9,
        "conditional_architecture": {
            "W0": "dimensionless strict informational operator and certified consequences",
            "CA": "energy/entropy/clock/length conversion resources",
            "SA": "state, orientation, gain, projective/spin and memory realization",
            "RA": "refinement and fine-data code",
            "OA": "preparation, channel, instrument, calibration, record and external custody",
        },
        "theorem": (
            "ST70-ST80 strengthen reproducibility and produce conditional source packages, but every proposed selector, "
            "refinement, thermodynamic, nonlinear, projective and custody mechanism still contains at least one added "
            "typed object. Removal countermodels therefore retain all nine source groups relative to the declared targets."
        ),
        "status": "proven_no_strict_source_group_reduction_after_ST70_ST80",
        "boundary": "The package architecture is a transparent conditional completion, not strict ToE closure or an absolute theorem over all future object classes.",
    }


def make_figures(results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    replay = results["ST70"]
    intervals = [replay["left_polynomial_interval"], replay["right_polynomial_interval"]]
    axes[0].bar(["left", "right"], [float(sum(map(Fraction, interval)) / 2) for interval in intervals])
    axes[0].axhline(0, color="black", lw=0.8); axes[0].set(title="exact-rational root signs", ylabel="p(s)")
    rows = results["ST71"]["rows"]
    axes[1].plot([row["odd_sector_curvature"] for row in rows], [row["gradient_flow_gain"] for row in rows], "o-")
    axes[1].axhline(0, color="black", lw=0.8); axes[1].set(title="oriented response classification", xlabel="odd Hessian", ylabel="gain")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st70_replay_st71_gain.png", dpi=190); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    st72 = results["ST72"]
    axes[0].bar(["rank delta", "rank augmented"], [st72["coboundary_matrix_rank"], st72["coboundary_augmented_rank"]])
    axes[0].set(title="nontrivial C24 cocycle")
    trajectory = results["ST73"]["online_trajectory"]
    axes[1].plot([row["time"] for row in trajectory], [row["supplied_q"] for row in trajectory], label="supplied")
    axes[1].plot([row["time"] for row in trajectory], [row["selected_q"] for row in trajectory], label="online selected")
    axes[1].set(title="fine-data refinement follows record", xlabel="step", ylabel="q"); axes[1].legend()
    fig.tight_layout(); fig.savefig(FIG_DIR / "st72_spin_st73_refinement.png", dpi=190); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    st74 = results["ST74"]
    ax.bar(["necessary", "sufficient", "ST51"], [st74["necessary_mean_shots_per_preparation_ceiling"], st74["distribution_free_sufficient_mean_shots_per_preparation_ceiling"], 1200])
    ax.set_yscale("log"); ax.set(title="nuisance-aware synthetic count bracket", ylabel="shots per preparation")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st74_robust_counts.png", dpi=190); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    st75 = results["ST75"]
    axes[0].bar(range(len(st75["spectral_block_multiplicities"])), st75["spectral_block_multiplicities"])
    axes[0].set(title="fixed equal-energy blocks", xlabel="block", ylabel="multiplicity")
    axes[1].bar([row["channel"] for row in st75["example_channels"]], [row["maximum_intertwining_residual"] for row in st75["example_channels"]])
    axes[1].set_yscale("log"); axes[1].set(title="CPTP examples", ylabel="intertwining residual")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st75_cptp_cone.png", dpi=190); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.2))
    resources = results["ST76"]["resources"]
    ax.barh(range(len(resources)), [1] * len(resources))
    ax.set_yticks(range(len(resources)), [row["resource"] for row in resources]); ax.set_xticks([])
    ax.set(title="independently necessary conditional thermodynamic resources")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st76_resource_minimality.png", dpi=190); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    branch = results["ST77"]["stored_branch_rows"]
    scatter = ax.scatter([row["kappa"] for row in branch], [row["IPR"] for row in branch], c=[row["tangent_kappa_component"] for row in branch], cmap="coolwarm", s=16)
    ax.axvline(results["ST77"]["fold_kappa"], color="black", ls=":", label="fold")
    ax.set(title="pseudo-arclength crosses the ST65 fold", xlabel="coupling kappa", ylabel="IPR"); ax.legend(); fig.colorbar(scatter, ax=ax, label="tangent kappa component")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st77_pseudo_arclength.png", dpi=190); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    st78 = results["ST78"]
    axes[0].plot([row["coupling"] for row in st78["rows"]], [row["positive_radial_root_count"] for row in st78["rows"]], "o-")
    axes[0].axvline(st78["critical_coupling_mu_over_lambda1"], color="red", ls=":")
    axes[0].set(title="strict stiffness quenches branches", xlabel="coupling", ylabel="positive radial roots")
    st79 = results["ST79"]
    axes[1].plot([row["time"] for row in st79["rows"]], [abs(row["strict_functional_calculus_chirality"]) for row in st79["rows"]], "o-")
    axes[1].set_yscale("log"); axes[1].set(title="strict Bargmann candidate vanishes", xlabel="time", ylabel="absolute chirality")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st78_backreaction_st79_chirality.png", dpi=190); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    st80 = results["ST80"]
    axes[0].bar(["valid", "tamper rejected", "duplicate rejected"], [int(st80["valid_case_accepted"]), int(st80["tampered_event_rejected"]), int(st80["duplicate_public_key_rejected"])])
    axes[0].set_ylim(0, 1.15); axes[0].tick_params(axis="x", rotation=15); axes[0].set(title="signed validator cases")
    st81 = results["ST81"]
    axes[1].bar(["before", "after"], [st81["strict_source_group_count_before"], st81["strict_source_group_count_after"]])
    axes[1].set(title="strict source groups remain", ylabel="count")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st80_validator_st81_sources.png", dpi=190); plt.close(fig)


def write_summary(results: dict) -> None:
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["program", "object", "status"])
        for number in range(70, 82):
            item = results[f"ST{number}"]
            writer.writerow([item["program"], item["object"], item["status"]])


def main() -> None:
    _, a, _ = strict_operator()
    results: dict[str, Any] = {
        "metadata": {
            "programs": "ST70-ST81",
            "seed": SEED,
            "date": "2026-08-11",
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sympy": sp.__version__,
            "cryptography": cryptography_version,
        }
    }
    results["ST70"] = st70_exact_rational_replay()
    results["ST71"] = st71_time_oriented_response_classification(a)
    results["ST72"] = st72_spin_central_extension_obstruction(a)
    results["ST73"] = st73_nonstationary_fine_data_refinement(a)
    results["ST74"] = st74_nuisance_robust_count_design(a)
    results["ST75"] = st75_all_cptp_intertwiners(a)
    results["ST76"] = st76_thermodynamic_resource_minimality()
    results["ST77"] = st77_pseudo_arclength_fold(a)
    results["ST78"] = st78_strict_doublet_backreaction(a)
    results["ST79"] = st79_orientation_odd_bargmann_obstruction(a)
    results["ST80"] = st80_signed_custody_validator(results["ST74"])
    results["ST81"] = st81_architecture_reconciliation()
    results["recommended_next_programs"] = [
        {"id": "ST82", "priority": 1, "study": "independently regenerate every ST58 transcendental and matrix enclosure with Arb/Sage or a proof assistant"},
        {"id": "ST83", "priority": 2, "study": "interval-certify the ST77 simple fold and local return branch in the seven-variable reflection sector"},
        {"id": "ST84", "priority": 3, "study": "seek a strict temporal-order source in admissible memory realizations or prove a passivity/detailed-balance obstruction"},
        {"id": "ST85", "priority": 4, "study": "classify state- or boundary-sourced lifts to the C24 central extension without silent selector insertion"},
        {"id": "ST86", "priority": 5, "study": "derive a fine-record statistic from an explicit local compression/RG map and test refinement universality"},
        {"id": "ST87", "priority": 6, "study": "replace the conservative ST74 concentration bound by a robust likelihood convex dual and exact nuisance optimization"},
        {"id": "ST88", "priority": 7, "study": "classify extreme points and reversible faces of the ST75 CPTP intertwiner cone"},
        {"id": "ST89", "priority": 8, "study": "formalize ST76 conversion-resource independence and distinguish thermal operations from Gibbs-preserving maps"},
        {"id": "ST90", "priority": 9, "study": "continue the ST77 branch beyond its first return and run a symmetry-classified disconnected-branch search"},
        {"id": "ST91", "priority": 10, "study": "analyze dynamic strict-doublet/order-parameter backreaction and bifurcation stability across the ST78 threshold"},
        {"id": "ST92", "priority": 11, "study": "search one genuinely strict-sourced reflection-odd object outside scalar functional calculus"},
        {"id": "ST93", "priority": 12, "study": "export a machine-readable W0+CA+SA+RA+OA conditional theory schema with claim-level dependency proofs"},
    ]
    results["epistemic_boundary"] = (
        "ST70-ST81 produce exact finite results, a partial independent replay, conditional constructions, numerical continuation, "
        "and synthetic validation. They do not supply QW-2191 closure, a strict time/gain/spin/refinement/dimensional source, "
        "laboratory evidence, legacy-to-strict completion or legacy role transfer, Standard Model, gravity, L_total, or ToE closure."
    )
    make_figures(results)
    write_summary(results)
    RESULTS.write_text(json.dumps(native(results), indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"results": RESULTS.name, "programs": 12, "figures": 8}, indent=2))


if __name__ == "__main__":
    main()
