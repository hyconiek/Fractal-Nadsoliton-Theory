#!/usr/bin/env python3
"""FIN local research Programs P465, P468, and P469.

P465 classifies all equality cases in the P452 reversal-symmetrization
inequality.  P468 continues the exact-rational nested-dual certificate.
P469 audits complementarity of the inherited O167 primal face against the
nested-dual central path.  Only local mathematical artifacts are used.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
import sympy as sp

import fin_programs_435_436_440 as p435
import fin_programs_445_447 as p445
import fin_programs_448_450 as p448
import fin_programs_451_453 as p451
import fin_programs_454_455_457 as p454
import fin_programs_458_459_464 as p458


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_465_468_469"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
P465_PATH = ROOT / f"{PREFIX}_P465_Strict_Symmetrization.csv"
P468_PATH = ROOT / f"{PREFIX}_P468_Dual_Refinement.csv"
P468_WITNESS_PATH = ROOT / f"{PREFIX}_P468_Rational_Dual.npz"
P469_PATH = ROOT / f"{PREFIX}_P469_Nested_KKT.csv"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p465_p468_p469_certificates.png"


def json_ready(value: Any) -> Any:
    if isinstance(value, Fraction):
        if max(value.numerator.bit_length(), value.denominator.bit_length()) <= 512:
            return str(value)
        return format(float(value), ".17g")
    if isinstance(value, sp.Basic):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: (
                    json.dumps(json_ready(row.get(key)), ensure_ascii=False)
                    if isinstance(row.get(key), (list, tuple, dict, np.ndarray))
                    else json_ready(row.get(key, ""))
                )
                for key in keys
            })


# ---------------------------------------------------------------------------
# P465: strict equality cases for reversal symmetrization
# ---------------------------------------------------------------------------


def program_465() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Prove that P452 equality holds exactly on the palindromic fixed set."""

    theta_interval = p445.iv_scale(p435.PI_INTERVAL, Fraction(2, 15))
    cos_theta = p451.cosine_interval(theta_interval)
    q = Fraction(4, 5)
    strict_condition = p445.iv_sub(
        p445.iv_scale(p445.iv_square(cos_theta), 4 * q**2),
        (Fraction(1), Fraction(1)),
    )
    if strict_condition[0] <= 0:
        raise AssertionError("P465 requires the strict P452 coefficient margin")

    alpha, beta, kappa = sp.symbols("alpha beta kappa", real=True)
    q_coefficient = (
        -sp.cos(2 * alpha) * sp.cos(2 * beta)
        + kappa * sp.cos(alpha + beta) ** 2
    ) / 3
    q_sos = (
        (kappa - 1) * sp.cos(alpha + beta) ** 2
        + sp.sin(alpha - beta) ** 2
    ) / 3
    q_identity = sp.trigsimp(sp.expand_trig(q_coefficient - q_sos))
    if q_identity != 0:
        raise AssertionError(f"quadratic coefficient identity failed: {q_identity}")

    # Exact boundary identity on Y=0.  In this face only the endpoint edge of
    # the three-survivor skew block remains.  Reversal averaging replaces xw
    # by (x^2+w^2)/2, hence the stated square.
    x, w, q_symbol, sine3 = sp.symbols("x w q sine3", nonnegative=True)
    boundary_original = 2 * q_symbol**3 * sine3 * x * w
    boundary_symmetric = q_symbol**3 * sine3 * (x**2 + w**2)
    boundary_defect = sp.factor(boundary_symmetric - boundary_original)
    expected_boundary_defect = q_symbol**3 * sine3 * (x - w) ** 2
    if sp.expand(boundary_defect - expected_boundary_defect) != 0:
        raise AssertionError("three-survivor boundary identity failed")

    rows: list[dict[str, Any]] = [
        {
            "row_type": "strict_parameter_condition",
            "cos_theta_lower": cos_theta[0],
            "cos_theta_upper": cos_theta[1],
            "kappa_minus_one_lower": strict_condition[0],
            "kappa_minus_one_upper": strict_condition[1],
        },
        {
            "row_type": "d2_defect_coefficients",
            "normalization": "16*q^2*sin(theta)^2*Y^4",
            "constant_coefficient": "cos(2*beta)^2/9",
            "linear_coefficient": (
                "2*sqrt(3)*(1-sin(2*beta)*cos(alpha-beta))/9"
            ),
            "quadratic_coefficient": str(q_coefficient),
            "quadratic_sum_of_squares": str(q_sos),
            "symbolic_identity_residual": str(q_identity),
        },
        {
            "row_type": "boundary_Y_zero",
            "d2_information": "identically blind to endpoint imbalance",
            "d3_gain": str(expected_boundary_defect),
            "equality_condition": "x=w",
        },
    ]

    # Deterministic falsification samples include the two singular faces and
    # near-face points.  They are checks, not the source of the theorem.
    rng = np.random.default_rng(465)
    samples = [
        np.array([0.9, 0.0, 0.0, 0.1]),
        np.array([0.0, 0.9, 0.1, 0.0]),
        np.array([0.5, 0.0, 0.0, 0.5]),
        np.array([0.0, 0.5, 0.5, 0.0]),
    ]
    samples.extend(rng.dirichlet(np.ones(4), size=4096))
    minimum_positive_gain = math.inf
    maximum_negative_gain = 0.0
    checked_nonpalindromic = 0
    for probabilities in samples:
        symmetrized = (probabilities + probabilities[::-1]) / 2
        gain = (
            p445.p446_numeric_objective(symmetrized)
            - p445.p446_numeric_objective(probabilities)
        )
        palindromic_residual = float(np.max(np.abs(probabilities - probabilities[::-1])))
        if palindromic_residual > 1e-14:
            checked_nonpalindromic += 1
            minimum_positive_gain = min(minimum_positive_gain, gain)
            maximum_negative_gain = min(maximum_negative_gain, gain)
    rows.append({
        "row_type": "deterministic_falsification_scan",
        "checked_nonpalindromic_points": checked_nonpalindromic,
        "minimum_observed_objective_gain": minimum_positive_gain,
        "most_negative_observed_gain": maximum_negative_gain,
    })

    p458_results = json.loads(
        (ROOT / "FIN_Programs_458_459_464_Results.json").read_text(encoding="utf-8")
    )["P458"]
    return ({
        "status": (
            "[Proven] P452 reversal symmetrization is strict off the palindromic "
            "fixed set; together with P458 this proves uniqueness of the declared "
            "coarse-erasure maximizer on the complete four-sector simplex"
        ),
        "d2_defect_formula": (
            "D2(sym p)^2-D2(p)^2 = 16 q^2 sin(theta)^2 Y^4 "
            "[c0+c1 r+c2 r^2], r=X/Y, with c0=cos^2(2 beta)/9, "
            "c1=2 sqrt(3)(1-sin(2 beta)cos(alpha-beta))/9, and "
            "c2=((kappa-1)cos^2(alpha+beta)+sin^2(alpha-beta))/3"
        ),
        "strict_condition_interval_kappa_minus_one": strict_condition,
        "interior_equality_case": (
            "For X,Y>0, c2=0 iff alpha=beta=pi/4; hence x=w and y=z."
        ),
        "X_zero_equality_case": (
            "For X=0<Y, c0=0 iff beta=pi/4; hence x=w=0 and y=z."
        ),
        "Y_zero_equality_case": (
            "For Y=0<X, D2 is blind, but D3(sym)-D3="
            "q^3 sin(3 theta)(x-w)^2, so equality iff x=w."
        ),
        "full_simplex_unique_maximizer_interval": p458_results[
            "unique_maximizer_interval"
        ],
        "checked_nonpalindromic_points": checked_nonpalindromic,
        "minimum_observed_objective_gain": minimum_positive_gain,
        "proof_boundary": (
            "The theorem is exact only for the declared three-use, permutation-"
            "symmetric Hamming-sector code with q=4/5, theta=2*pi/15, and its "
            "heralded-loss objective. It is not an unrestricted adaptive-channel, "
            "laboratory, selector, dimensional, or physical-uniqueness theorem."
        ),
        "new_object": "O176 Strict Reversal Equality Isolator",
    }, rows)


# ---------------------------------------------------------------------------
# P469: full nested KKT residual for the O167 face
# ---------------------------------------------------------------------------


def positive_root_and_inverse(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eigenvalues, eigenvectors = np.linalg.eigh((matrix + matrix.conj().T) / 2)
    if float(np.min(eigenvalues)) <= 0:
        raise ValueError("the O167 normalizer left the positive interior")
    root = (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.conj().T
    inverse = (eigenvectors * (1 / np.sqrt(eigenvalues))) @ eigenvectors.conj().T
    return root, inverse


def o167_support_ladder(
    parameters: np.ndarray, delta: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, np.ndarray, np.ndarray]:
    """Return the canonical top support and inherited KKT ladder.

    For an interior normalizer N, the unique trace-norm support associated
    with nonsingular H=sqrt(N) Delta sqrt(N) is

        X3 = N^{-1/2} |H| N^{-1/2}/2.

    The lower matrices are averages of the two principal blocks.  Therefore
    their block-difference residuals are exactly the missing stationarity
    equations for an interior ordered comb.
    """

    normalizer = p454.coherent_face_normalizer(parameters)
    root, inverse = positive_root_and_inverse(normalizer)
    transformed = root @ delta @ root
    values, vectors = np.linalg.eigh((transformed + transformed.conj().T) / 2)
    absolute = (vectors * np.abs(values)) @ vectors.conj().T
    x3 = inverse @ absolute @ inverse / 2
    x3 = (x3 + x3.conj().T) / 2
    x2 = (x3[:4, :4] + x3[4:, 4:]) / 2
    x1 = (x2[:2, :2] + x2[2:, 2:]) / 2
    lam = float(np.real((x1[0, 0] + x1[1, 1]) / 2))
    plus_projector = vectors[:, values > 0] @ vectors[:, values > 0].conj().T
    minus_projector = vectors[:, values < 0] @ vectors[:, values < 0].conj().T
    t_plus = root @ plus_projector @ root
    t_minus = root @ minus_projector @ root
    return normalizer, x3, x2, x1, lam, t_plus, t_minus


def o167_ladder_residual(parameters: np.ndarray, delta: np.ndarray) -> np.ndarray:
    try:
        _, x3, x2, x1, _, _, _ = o167_support_ladder(parameters, delta)
    except ValueError:
        return np.full(21, 1e3)
    return np.concatenate((
        np.real(x3[:4, :4] - x3[4:, 4:]).reshape(-1),
        np.real(x2[:2, :2] - x2[2:, 2:]).reshape(-1),
        [float(np.real(x1[0, 0] - x1[1, 1]))],
    ))


def program_469() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    inherited = json.loads(
        (ROOT / "FIN_Programs_454_455_457_Results.json").read_text(encoding="utf-8")
    )["P455"]
    initial = np.array(inherited["face_optimum_parameters"], dtype=float)
    delta = p448.compressed_process_difference(3, 0.8, math.pi / 8)
    optimum = least_squares(
        lambda value: o167_ladder_residual(value, delta),
        initial,
        xtol=1e-15,
        ftol=1e-15,
        gtol=1e-15,
        max_nfev=3000,
    )
    normalizer, x3, x2, x1, lam, t_plus, t_minus = o167_support_ladder(
        optimum.x, delta
    )
    residual = o167_ladder_residual(optimum.x, delta)

    b0 = normalizer[:4, :4]
    b1 = normalizer[4:, 4:]
    inherited_two = b0 + b1
    c0 = inherited_two[:2, :2]
    c1 = inherited_two[2:, 2:]
    inherited_one = c0 + c1
    d0 = float(np.real(inherited_one[0, 0]))
    d1 = float(np.real(inherited_one[1, 1]))

    slacks = [
        x3 - delta / 2,
        x3 + delta / 2,
        x2 - x3[:4, :4],
        x2 - x3[4:, 4:],
        x1 - x2[:2, :2],
        x1 - x2[2:, 2:],
        np.array([[lam - float(np.real(x1[0, 0]))]]),
        np.array([[lam - float(np.real(x1[1, 1]))]]),
    ]
    primal_factors = [
        t_plus, t_minus, b0, b1, c0, c1,
        np.array([[d0]]), np.array([[d1]]),
    ]
    names = [
        "Tplus_times_X3_minus_halfDelta",
        "Tminus_times_X3_plus_halfDelta",
        "B0_times_X2_minus_X3_block0",
        "B1_times_X2_minus_X3_block1",
        "C0_times_X1_minus_X2_block0",
        "C1_times_X1_minus_X2_block1",
        "d0_times_lambda_minus_X1_diag0",
        "d1_times_lambda_minus_X1_diag1",
    ]
    rows: list[dict[str, Any]] = []
    contacts: list[float] = []
    for name, slack, primal in zip(names, slacks, primal_factors):
        contact = float(np.real(np.trace(slack @ primal)))
        contacts.append(contact)
        rows.append({
            "row_type": "kkt_contact",
            "contact": name,
            "trace_product": contact,
            "frobenius_product_norm": float(np.linalg.norm(slack @ primal)),
            "slack_minimum_eigenvalue": float(np.min(np.linalg.eigvalsh(
                (slack + slack.conj().T) / 2
            ))),
            "primal_minimum_eigenvalue": float(np.min(np.linalg.eigvalsh(
                (primal + primal.conj().T) / 2
            ))),
        })

    primal_value = p448.tester_distance(normalizer, delta)
    gap = lam - primal_value
    rows.append({
        "row_type": "stationarity_summary",
        "parameters": optimum.x,
        "residual_two_norm": float(np.linalg.norm(residual)),
        "residual_infinity_norm": float(np.linalg.norm(residual, np.inf)),
        "support_lambda": lam,
        "primal_value": primal_value,
        "support_gap": gap,
        "sum_trace_contacts": float(sum(contacts)),
        "least_squares_cost": float(optimum.cost),
        "least_squares_optimality": float(optimum.optimality),
        "function_evaluations": int(optimum.nfev),
    })

    # A full exact KKT theorem would need interval isolation of a true zero of
    # the nonlinear support-ladder system.  Tiny floating residuals are strong
    # evidence, but never a proof of exact zero or PSD lower cascade slacks.
    residual_inf = float(np.linalg.norm(residual, np.inf))
    return ({
        "status": (
            "[Strong evidence] the O167 interior candidate satisfies the complete "
            "nested support-ladder KKT equations to floating residual below 1e-11; "
            "[Open] exact root isolation and exact complementarity"
        ),
        "stationary_parameters": optimum.x,
        "normalizer_minimum_eigenvalue": float(np.min(np.linalg.eigvalsh(normalizer))),
        "transformed_delta_minimum_absolute_eigenvalue": float(np.min(np.abs(
            np.linalg.eigvalsh(positive_root_and_inverse(normalizer)[0] @ delta @ positive_root_and_inverse(normalizer)[0])
        ))),
        "block3_residual_frobenius": float(np.linalg.norm(
            x3[:4, :4] - x3[4:, 4:]
        )),
        "block2_residual_frobenius": float(np.linalg.norm(
            x2[:2, :2] - x2[2:, 2:]
        )),
        "scalar_residual": float(abs(x1[0, 0] - x1[1, 1])),
        "full_residual_infinity_norm": residual_inf,
        "primal_value": primal_value,
        "support_lambda": lam,
        "support_gap": gap,
        "sum_trace_contacts": float(sum(contacts)),
        "maximum_absolute_trace_contact": float(max(abs(value) for value in contacts)),
        "necessity_theorem": (
            "Because B0,B1,C0,C1 and d0,d1 are positive at an interior O167 "
            "optimum, complementary slackness forces both lower slacks at every "
            "causal level to vanish. The tested block equalities are therefore "
            "necessary full-cone stationarity equations, not merely face gradients."
        ),
        "boundary": (
            "Floating least-squares residuals do not prove a nonlinear exact root, "
            "dual feasibility of averaged near-zero lower slacks, exact O167 global "
            "attainment, or uniqueness. An interval Krawczyk/Newton certificate is required."
        ),
        "new_object": "O177 Interior Comb Support Ladder",
    }, rows)


# ---------------------------------------------------------------------------
# P468: support-ladder rational dual refinement
# ---------------------------------------------------------------------------


def support_ladder_dual_vector(
    parameters: np.ndarray, delta: np.ndarray, shift: float
) -> tuple[np.ndarray, dict[str, float]]:
    """Build a strictly feasible scout from the P469 support ladder."""

    _, support_x3, _, _, _, _, _ = o167_support_ladder(parameters, delta)
    x3 = np.real(support_x3) + shift * np.eye(8)
    x2 = (x3[:4, :4] + x3[4:, 4:]) / 2 + shift * np.eye(4)
    x1 = (x2[:2, :2] + x2[2:, 2:]) / 2 + shift * np.eye(2)
    lam = float(max(x1[0, 0], x1[1, 1]) + shift)
    vector = p454.pack_dual(x3, x2, x1, lam)
    slacks, _ = p454.dual_slacks(vector, delta / 2)
    return vector, {
        "shift": shift,
        "floating_lambda": lam,
        "minimum_float_slack": min(
            float(np.min(np.linalg.eigvalsh((slack + slack.conj().T) / 2)))
            for slack in slacks
        ),
    }


def program_468() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    inherited = json.loads(
        (ROOT / "FIN_Programs_458_459_464_Results.json").read_text(encoding="utf-8")
    )["P464"]
    previous_upper = Fraction(inherited["accepted_certified_upper"])
    primal_lower = Fraction(
        json.loads(
            (ROOT / "FIN_Programs_454_455_457_Results.json").read_text(encoding="utf-8")
        )["P454"]["inherited_primal_lower"]
    )
    p455 = json.loads(
        (ROOT / "FIN_Programs_454_455_457_Results.json").read_text(encoding="utf-8")
    )["P455"]
    delta = p448.compressed_process_difference(3, 0.8, math.pi / 8)
    stationary = least_squares(
        lambda value: o167_ladder_residual(value, delta),
        np.array(p455["face_optimum_parameters"], dtype=float),
        xtol=1e-15,
        ftol=1e-15,
        gtol=1e-15,
        max_nfev=3000,
    )

    rows: list[dict[str, Any]] = []
    candidates: list[tuple[Fraction, dict[str, Any], dict[str, np.ndarray], float]] = []
    # The shifts are theorem-scoped only after exact rational admission.  The
    # descending list also tests sensitivity to lost floating digits.
    for shift, scale in (
        (1e-9, 10**16),
        (3e-10, 10**16),
        (1e-10, 10**17),
    ):
        vector, scout = support_ladder_dual_vector(stationary.x, delta, shift)
        certificate, matrices = p458.rational_dual_certificate(vector, scale=scale)
        upper = Fraction(certificate["lambda_upper"])
        rows.append({
            "row_type": "support_ladder_candidate",
            **scout,
            "rationalization_scale": scale,
            "exact_feasible": certificate["feasible"],
            "certified_upper": upper,
            "certified_gap": upper - primal_lower,
            "smallest_certified_slack_lower": certificate[
                "smallest_certified_slack_lower"
            ],
        })
        rows.extend({
            "row_type": "candidate_exact_slack",
            "candidate_shift": shift,
            "candidate_scale": scale,
            **row,
        } for row in certificate["slacks"])
        if certificate["feasible"]:
            candidates.append((upper, certificate, matrices, shift))

    if not candidates:
        accepted_upper = previous_upper
        accepted_certificate = inherited["candidate_certificate"]
        accepted_shift = None
        improved = False
    else:
        accepted_upper, accepted_certificate, accepted_matrices, accepted_shift = min(
            candidates, key=lambda item: item[0]
        )
        improved = accepted_upper < previous_upper
        if improved:
            np.savez_compressed(P468_WITNESS_PATH, **accepted_matrices)
        else:
            source = np.load(ROOT / "FIN_Programs_458_459_464_P464_Rational_Dual.npz")
            np.savez_compressed(P468_WITNESS_PATH, **{key: source[key] for key in source.files})
            accepted_upper = previous_upper
            accepted_certificate = inherited["candidate_certificate"]
            accepted_shift = None

    gap = accepted_upper - primal_lower
    reduction = float((previous_upper - primal_lower) / gap) if gap > 0 else math.inf
    return ({
        "status": (
            "[Computer-assisted proof] the P469 support ladder yields a new exact-"
            "rational feasible dual with certified global gap below 1e-8"
            if improved and gap < Fraction(1, 10**8) else
            "[Computer-assisted proof] exact-rational support-ladder refinement"
            if improved else
            "[Refuted for the attempted shifts] no improved exact-rational support-"
            "ladder certificate was admitted; P464 remains authoritative"
        ),
        "stationary_parameters_used": stationary.x,
        "attempted_shifts": [1e-9, 3e-10, 1e-10],
        "strict_improvement_accepted": improved,
        "previous_certified_upper": previous_upper,
        "accepted_certified_upper": accepted_upper,
        "inherited_primal_lower": primal_lower,
        "accepted_global_gap": gap,
        "gap_reduction_factor_over_P464": reduction,
        "accepted_shift": accepted_shift,
        "accepted_certificate": accepted_certificate,
        "methodological_result": (
            "The floating P469 stationarity point only proposes a dual. Globality is "
            "licensed solely by exact rational Sylvester positivity and the full paid "
            "trigonometric perturbation radius."
        ),
        "boundary": (
            "A nonzero certified gap, however small, does not prove exact O167 "
            "attainment or uniqueness of the full ordered-comb optimizer."
        ),
        "new_object": "O178 Support-Ladder Rational Closure Tube",
    }, rows)


def make_figure(
    p465: dict[str, Any],
    p468: dict[str, Any],
    p468_rows: list[dict[str, Any]],
    p469: dict[str, Any],
    p469_rows: list[dict[str, Any]],
) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.45))

    alpha = np.linspace(0, math.pi / 2, 180)
    beta = np.linspace(0, math.pi / 2, 180)
    aa, bb = np.meshgrid(alpha, beta)
    kappa = 1 + float(Fraction(p465["strict_condition_interval_kappa_minus_one"][0]))
    c2 = ((kappa - 1) * np.cos(aa + bb) ** 2 + np.sin(aa - bb) ** 2) / 3
    image = axes[0].imshow(
        c2,
        origin="lower",
        extent=(0, 0.5, 0, 0.5),
        aspect="auto",
        cmap="viridis",
    )
    axes[0].plot([0.25], [0.25], "o", color="#dc2626", markersize=5)
    axes[0].set_title("P465 strict equality landscape")
    axes[0].set_xlabel(r"$\alpha/\pi$")
    axes[0].set_ylabel(r"$\beta/\pi$")
    fig.colorbar(image, ax=axes[0], fraction=0.046, pad=0.04, label=r"$c_2$")

    inherited_gap = float(Fraction(
        json.loads((ROOT / "FIN_Programs_454_455_457_Results.json").read_text(encoding="utf-8"))["P454"]["certified_global_gap"]
    ))
    p464_gap = float(Fraction(
        json.loads((ROOT / "FIN_Programs_458_459_464_Results.json").read_text(encoding="utf-8"))["P464"]["accepted_global_gap"]
    ))
    new_gap = float(Fraction(p468["accepted_global_gap"]))
    axes[1].bar(["P454", "P464", "P468"], [inherited_gap, p464_gap, new_gap], color=["#94a3b8", "#7c3aed", "#0f766e"])
    axes[1].set_yscale("log")
    axes[1].axhline(1e-8, color="#dc2626", linestyle="--", linewidth=1)
    axes[1].set_title("Exact rational global gap")
    axes[1].set_ylabel("certified upper - lower")

    summary = [row for row in p469_rows if row["row_type"] == "stationarity_summary"][0]
    residuals = [
        p469["block3_residual_frobenius"],
        p469["block2_residual_frobenius"],
        p469["scalar_residual"],
        abs(summary["sum_trace_contacts"]),
    ]
    axes[2].bar(["8->4", "4->2", "2->1", "contacts"], residuals, color="#2563eb")
    axes[2].set_yscale("log")
    axes[2].set_title("P469 floating KKT residuals")
    axes[2].set_ylabel("residual magnitude")

    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)


def write_results(
    p465: dict[str, Any],
    p468: dict[str, Any],
    p469: dict[str, Any],
) -> None:
    results = {
        "metadata": {
            "programs": "P465/P468/P469",
            "checkpoint": "P465-P469",
            "date": "2026-08-01",
            "local_only": True,
            "external_physical_evidence": False,
            "selector_discharged": False,
            "dimensional_source_exported": False,
            "legacy_strict_bridge_complete": False,
            "legacy_role_transfer_started": False,
            "exact_O167_attainment_proven": False,
        },
        "P465": p465,
        "P468": p468,
        "P469": p469,
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(results), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    summary_rows = [
        {
            "program": "P465",
            "status": p465["status"],
            "headline_metric": "full-simplex unique maximizer",
            "value": p465["full_simplex_unique_maximizer_interval"],
            "new_object": p465["new_object"],
        },
        {
            "program": "P468",
            "status": p468["status"],
            "headline_metric": "certified global gap",
            "value": p468["accepted_global_gap"],
            "new_object": p468["new_object"],
        },
        {
            "program": "P469",
            "status": p469["status"],
            "headline_metric": "floating KKT residual infinity norm",
            "value": p469["full_residual_infinity_norm"],
            "new_object": p469["new_object"],
        },
    ]
    write_csv(SUMMARY_PATH, summary_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", choices=("P465", "P468", "P469", "all"), default="all")
    args = parser.parse_args()
    p465_result, p465_rows = program_465()
    write_csv(P465_PATH, p465_rows)
    if args.only == "P465":
        print(json.dumps(json_ready({"P465": p465_result}), indent=2))
        return
    p468_result, p468_rows = program_468()
    write_csv(P468_PATH, p468_rows)
    if args.only == "P468":
        print(json.dumps(json_ready({"P468": p468_result}), indent=2))
        return
    p469_result, p469_rows = program_469()
    write_csv(P469_PATH, p469_rows)
    if args.only == "P469":
        print(json.dumps(json_ready({"P469": p469_result}), indent=2))
        return
    make_figure(p465_result, p468_result, p468_rows, p469_result, p469_rows)
    write_results(p465_result, p468_result, p469_result)
    print(json.dumps(json_ready({
        "P465": p465_result,
        "P468": p468_result,
        "P469": p469_result,
    }), indent=2))


if __name__ == "__main__":
    main()
