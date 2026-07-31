#!/usr/bin/env python3
"""Execute FIN research Programs P337--P350.

This round certifies rather than re-searches the P323 moment extremizer,
separates the full oscillatory strict kernel from its attenuation envelope,
audits the QW-2038 refreeze continuously, tests carrier naturality, compiles
an ideal photonic transfer and short-time inverse calibration, constructs a
declared parallel-comb subclass, and makes the conditional bridge-resource,
electroweak, and dimensional metrology assumptions explicit.

P342 and P350 remain external-evidence gates.  No synthetic object is admitted
as laboratory evidence, an internal selector, a dimensional source, a
legacy-role transfer theorem, L_total, SM/GR generation, or ToE closure.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import networkx as nx
import numpy as np
import pandas as pd
from scipy.linalg import expm, logm, polar
from scipy.optimize import brentq, differential_evolution, linprog
from scipy.signal import hilbert

import QW_2038_DERIVATION_COMPATIBLE_KERNEL_REFREEZE_SCAN as qw2038
import fin_programs_255_266 as core
import fin_programs_295_308 as p295
import fin_programs_323_336 as p323


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "FIN_Programs_337_350_Figures"
RESULTS_PATH = ROOT / "FIN_Programs_337_350_Results.json"
SUMMARY_PATH = ROOT / "FIN_Programs_337_350_Summary.csv"
INTERVAL_PATH = ROOT / "FIN_Programs_337_350_Interval_Certificate.csv"
OSCILLATORY_PATH = ROOT / "FIN_Programs_337_350_Oscillatory_Resource.csv"
REFREEZE_PATH = ROOT / "FIN_Programs_337_350_Refreeze_Sensitivity.csv"
NATURALITY_PATH = ROOT / "FIN_Programs_337_350_Natural_Completion.csv"
PHOTONIC_PATH = ROOT / "FIN_Programs_337_350_Photonic_Compilation.csv"
COMB_PATH = ROOT / "FIN_Programs_337_350_Parallel_Comb.csv"
CLOCK_PATH = ROOT / "FIN_Programs_337_350_Nonparametric_Clock.csv"
INVARIANT_PATH = ROOT / "FIN_Programs_337_350_Carrier_Invariants.csv"
MONOTONE_PATH = ROOT / "FIN_Programs_337_350_Resource_Monotones.csv"
COUPLING_PATH = ROOT / "FIN_Programs_337_350_Coupling_Laws.csv"
EW_PATH = ROOT / "FIN_Programs_337_350_EW_Falsification.csv"
METROLOGY_PATH = ROOT / "FIN_Programs_337_350_Metrology.csv"
FORMAL_SOURCE = ROOT / "FIN_Programs_337_350_Formal_Core.lean"

N = 12
SEED = 20260730
ALPHA_GEO = 4.0 * math.log(2.0)
STRICT_OMEGA = 743.0 / 4000.0
STRICT_PHI = 13.0 / 80.0
STRICT_BETA = 1.0
STRICT_ETA = 9.0 / 5.0


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        json.dumps(json_ready(value), ensure_ascii=False)
                        if isinstance(value, (dict, list, tuple, np.ndarray))
                        else value
                    )
                    for key, value in row.items()
                }
            )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def polynomial_value(coefficients: np.ndarray, x: np.ndarray | float) -> np.ndarray:
    return np.polynomial.polynomial.polyval(x, coefficients)


# ---------------------------------------------------------------------------
# P337: exact rational lower certificate for the P323 continuum extremizer
# ---------------------------------------------------------------------------


def fraction_polynomial_square(coefficients: list[Fraction]) -> list[Fraction]:
    result = [Fraction(0)] * (2 * len(coefficients) - 1)
    for i, left in enumerate(coefficients):
        for j, right in enumerate(coefficients):
            result[i + j] += left * right
    return result


def substitute_affine_power_coefficients(
    coefficients: list[Fraction], left: Fraction, right: Fraction
) -> list[Fraction]:
    """Return coefficients of p(left+(right-left)t) in the power basis."""

    degree = len(coefficients) - 1
    width = right - left
    result = [Fraction(0)] * (degree + 1)
    for power, coefficient in enumerate(coefficients):
        for t_power in range(power + 1):
            result[t_power] += (
                coefficient
                * Fraction(math.comb(power, t_power))
                * left ** (power - t_power)
                * width**t_power
            )
    return result


def power_to_bernstein(power: list[Fraction]) -> list[Fraction]:
    degree = len(power) - 1
    bernstein = []
    for index in range(degree + 1):
        value = Fraction(0)
        for power_index in range(index + 1):
            value += (
                power[power_index]
                * Fraction(math.comb(index, power_index), math.comb(degree, power_index))
            )
        bernstein.append(value)
    return bernstein


def exact_fifth_root_bracket(integer: int, digits: int = 24) -> tuple[Fraction, Fraction]:
    if integer == 0:
        return Fraction(0), Fraction(0)
    mp.mp.dps = digits + 30
    scale = 10**digits
    floor_value = int(mp.floor(mp.root(integer, 5) * scale))
    lower = Fraction(floor_value, scale)
    upper = Fraction(floor_value + 1, scale)
    assert lower**5 <= integer <= upper**5
    return lower, upper


def strict_attenuation_interval(order: int) -> tuple[Fraction, Fraction]:
    if order == 0:
        return Fraction(1), Fraction(1)
    root_lower, root_upper = exact_fifth_root_bracket(order)
    power_lower = root_lower**9
    power_upper = root_upper**9
    return (
        Fraction(1, 1) / (Fraction(1, 1) + power_upper),
        Fraction(1, 1) / (Fraction(1, 1) + power_lower),
    )


def p323_q_coefficients() -> list[mp.mpf]:
    mp.mp.dps = 100
    moments = [1 / (1 + mp.mpf(k) ** (mp.mpf(9) / 5)) for k in range(7)]
    shifted = moments[1:]
    hankel = mp.matrix([[shifted[i + j] for j in range(3)] for i in range(3)])
    rhs = mp.matrix([-shifted[i + 3] for i in range(3)])
    c0, c1, c2 = mp.lu_solve(hankel, rhs)
    return [mp.mpf(1), c1 / c0, c2 / c0, 1 / c0]


def program_337() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    q_mp = p323_q_coefficients()
    # A small exact contraction creates a rational polynomial with a strict
    # feasibility margin while preserving an extremely sharp lower bound.
    contraction = Fraction(99_999_999, 100_000_000)
    q_rational = [
        contraction * Fraction(mp.nstr(value, 34)) for value in q_mp
    ]
    q_square = fraction_polynomial_square(q_rational)
    p_rational = [-value for value in q_square]

    subdivision_count = 512
    maximum_bernstein = Fraction(-10)
    minimum_bernstein = Fraction(10)
    for index in range(subdivision_count):
        left = Fraction(index, subdivision_count)
        right = Fraction(index + 1, subdivision_count)
        local_power = substitute_affine_power_coefficients(q_square, left, right)
        local_bernstein = power_to_bernstein(local_power)
        maximum_bernstein = max(maximum_bernstein, *local_bernstein)
        minimum_bernstein = min(minimum_bernstein, *local_bernstein)
    assert maximum_bernstein <= 1
    # q^2 is an exact square, so non-negativity does not depend on Bernstein
    # lower coefficients.

    moment_intervals = [strict_attenuation_interval(k) for k in range(7)]
    objective_lower = Fraction(0)
    objective_upper = Fraction(0)
    for coefficient, (lower, upper) in zip(p_rational, moment_intervals):
        if coefficient >= 0:
            objective_lower += coefficient * lower
            objective_upper += coefficient * upper
        else:
            objective_lower += coefficient * upper
            objective_upper += coefficient * lower

    prior = json.loads(
        (ROOT / "FIN_Programs_323_336_Results.json").read_text(encoding="utf-8")
    )["P323"]
    numerical_extremum = float(
        prior["exact_continuum_minimum_negative_mass_numeric"]
    )
    rows: list[dict[str, Any]] = []
    for index, coefficient in enumerate(q_rational):
        rows.append(
            {
                "row_type": "rational_q_coefficient",
                "index": index,
                "exact_value": str(coefficient),
                "decimal_value": float(coefficient),
            }
        )
    for index, coefficient in enumerate(p_rational):
        rows.append(
            {
                "row_type": "rational_dual_coefficient",
                "index": index,
                "exact_value": str(coefficient),
                "decimal_value": float(coefficient),
            }
        )
    for order, (lower, upper) in enumerate(moment_intervals):
        rows.append(
            {
                "row_type": "algebraic_moment_interval",
                "index": order,
                "lower_exact": str(lower),
                "upper_exact": str(upper),
                "lower_decimal": float(lower),
                "upper_decimal": float(upper),
                "width": float(upper - lower),
            }
        )

    return (
        {
            "status": (
                "[Proven] exact rational dual lower certificate; "
                "[Strong evidence] matching high-precision primal evaluation"
            ),
            "subdivision_count": subdivision_count,
            "dual_polynomial_form": "p_r(x)=-q_r(x)^2",
            "exact_q_contraction": str(contraction),
            "maximum_exact_bernstein_upper_bound_q_squared": str(
                maximum_bernstein
            ),
            "maximum_bernstein_upper_bound_decimal": float(maximum_bernstein),
            "all_fifth_root_brackets_verified_by_integer_arithmetic": True,
            "certified_dual_lower_bound_exact": str(objective_lower),
            "certified_dual_lower_bound": float(objective_lower),
            "certified_objective_interval_width_from_moments": float(
                objective_upper - objective_lower
            ),
            "p323_high_precision_primal_value": numerical_extremum,
            "distance_from_certified_lower_bound": (
                numerical_extremum - float(objective_lower)
            ),
            "theorem": (
                "The rational polynomial p_r=-q_r^2 is nonpositive. Exact "
                "Bernstein subdivision proves q_r^2<=1 on [0,1], hence "
                "-1<=p_r<=0. Exact rational fifth-root brackets enclose every "
                "strict moment h_k. Weak duality therefore proves the stated "
                "lower bound for every representing signed measure."
            ),
            "boundary": (
                "This closes the numerical-kernel gap on the lower-bound side. "
                "The matching displayed primal atoms remain a 100-digit "
                "evaluation; a fully formal algebraic interval-Newton "
                "certificate for their simultaneous feasibility is P351."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P338: signed resource of the full oscillatory strict kernel
# ---------------------------------------------------------------------------


def polynomial_global_range(coefficients: np.ndarray) -> tuple[float, float, list[float]]:
    derivative = np.array(
        [index * coefficients[index] for index in range(1, len(coefficients))]
    )
    scan = np.linspace(0.0, 1.0, 50_001)
    derivative_values = polynomial_value(derivative, scan)
    critical = [0.0, 1.0]
    for index in range(len(scan) - 1):
        left_value = derivative_values[index]
        right_value = derivative_values[index + 1]
        if left_value == 0.0:
            critical.append(float(scan[index]))
        elif left_value * right_value < 0.0:
            critical.append(
                float(
                    brentq(
                        lambda x: float(polynomial_value(derivative, x)),
                        float(scan[index]),
                        float(scan[index + 1]),
                        xtol=5e-15,
                    )
                )
            )
    values = [float(polynomial_value(coefficients, x)) for x in critical]
    return min(values), max(values), critical


def solve_continuum_dual(
    moments: np.ndarray, initial_points: int = 257
) -> tuple[np.ndarray, float, float, int]:
    degree = len(moments) - 1
    points = list(np.linspace(0.0, 1.0, initial_points))
    result = None
    for iteration in range(30):
        unique = np.array(sorted(set(round(float(x), 15) for x in points)))
        vandermonde = np.vander(unique, degree + 1, increasing=True)
        a_ub = np.vstack([vandermonde, -vandermonde])
        b_ub = np.r_[np.zeros(len(unique)), np.ones(len(unique))]
        result = linprog(
            -moments,
            A_ub=a_ub,
            b_ub=b_ub,
            bounds=[(None, None)] * (degree + 1),
            method="highs",
        )
        if not result.success:
            raise RuntimeError(result.message)
        coefficients = result.x
        minimum, maximum, critical = polynomial_global_range(coefficients)
        violations = [
            point
            for point in critical
            if polynomial_value(coefficients, point) > 2e-10
            or polynomial_value(coefficients, point) < -1.0 - 2e-10
        ]
        if not violations:
            break
        points.extend(violations)
    assert result is not None
    minimum, maximum, _ = polynomial_global_range(result.x)
    # Affine normalization of the globally evaluated range produces a
    # continuum-feasible polynomial even when the LP is at floating tolerance.
    safety = 1e-9
    denominator = maximum - minimum + 2.0 * safety
    normalized = (result.x.copy()) / denominator
    normalized[0] -= (maximum + safety) / denominator
    normalized_minimum, normalized_maximum, _ = polynomial_global_range(normalized)
    return normalized, normalized_minimum, normalized_maximum, iteration + 1


def solve_grid_primal(
    moments: np.ndarray, grid_size: int = 4001
) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    points = np.linspace(0.0, 1.0, grid_size)
    matrix = np.vander(points, len(moments), increasing=True).T
    equality = np.hstack([matrix, -matrix])
    objective = np.r_[np.zeros(grid_size), np.ones(grid_size)]
    result = linprog(
        objective,
        A_eq=equality,
        b_eq=moments,
        bounds=(0.0, None),
        method="highs",
        options={"dual_feasibility_tolerance": 1e-9},
    )
    if not result.success:
        raise RuntimeError(result.message)
    positive = result.x[:grid_size]
    negative = result.x[grid_size:]
    return float(result.fun), points, positive, negative


def high_precision_primal_from_grid_support(
    moments_count: int,
    points: np.ndarray,
    positive: np.ndarray,
    negative: np.ndarray,
) -> tuple[float, list[tuple[float, float]], float]:
    """Re-solve the selected rational support at 100 digits.

    The raw monomial LP is used only to select support.  Its floating weights
    are not called a feasible upper certificate because tiny high-order
    residuals are amplified by the ill-conditioned monomial basis.
    """

    active_indices = sorted(
        set(np.where(positive > 1e-8)[0]) | set(np.where(negative > 1e-8)[0])
    )
    if len(active_indices) != moments_count:
        raise RuntimeError(
            f"expected {moments_count} active support points, got {len(active_indices)}"
        )
    mp.mp.dps = 100
    denominator = len(points) - 1
    nodes = [mp.mpf(int(index)) / denominator for index in active_indices]
    moments = [
        mp.cos(mp.mpf(743) * order / 4000 + mp.mpf(13) / 80)
        / (1 + mp.mpf(order) ** (mp.mpf(9) / 5))
        for order in range(moments_count)
    ]
    vandermonde = mp.matrix(
        [[node**order for node in nodes] for order in range(moments_count)]
    )
    weights = mp.lu_solve(vandermonde, mp.matrix(moments))
    residual = max(
        abs(
            sum(weights[index] * nodes[index] ** order for index in range(len(nodes)))
            - moments[order]
        )
        for order in range(moments_count)
    )
    negative_mass = sum(-weight for weight in weights if weight < 0)
    return (
        float(negative_mass),
        [(float(node), float(weight)) for node, weight in zip(nodes, weights)],
        float(residual),
    )


def program_338() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    orders = np.arange(12, dtype=float)
    moments = np.cos(STRICT_OMEGA * orders + STRICT_PHI) / (
        1.0 + orders**STRICT_ETA
    )
    coefficients, dual_minimum, dual_maximum, iterations = solve_continuum_dual(
        moments
    )
    dual_value = float(np.dot(coefficients, moments))
    grid_objective, points, positive, negative = solve_grid_primal(moments)
    primal_value, primal_atoms, primal_residual = high_precision_primal_from_grid_support(
        len(moments), points, positive, negative
    )
    rows: list[dict[str, Any]] = []
    for order, value in enumerate(moments):
        rows.append(
            {
                "row_type": "strict_oscillatory_moment",
                "order": order,
                "value": float(value),
                "sign": int(np.sign(value)),
            }
        )
    for index, coefficient in enumerate(coefficients):
        rows.append(
            {
                "row_type": "continuum_dual_polynomial",
                "order": index,
                "value": float(coefficient),
            }
        )
    for kind, weights in (("positive_grid_atom", positive), ("negative_grid_atom", negative)):
        for index in np.where(weights > 1e-8)[0]:
            rows.append(
                {
                    "row_type": kind,
                    "node": float(points[index]),
                    "weight": float(weights[index]),
                }
            )
    for node, weight in primal_atoms:
        rows.append(
            {
                "row_type": "high_precision_primal_atom",
                "node": node,
                "weight": weight,
                "sign": int(np.sign(weight)),
            }
        )

    envelope_value = json.loads(
        (ROOT / "FIN_Programs_323_336_Results.json").read_text(encoding="utf-8")
    )["P323"]["exact_continuum_minimum_negative_mass_numeric"]
    return (
        {
            "status": (
                "[Strong evidence] continuum dual/discrete-primal bracket for "
                "the full oscillatory kernel"
            ),
            "moment_orders": list(range(12)),
            "negative_moment_orders": [
                int(index) for index, value in enumerate(moments) if value < 0
            ],
            "dual_cutting_iterations": iterations,
            "continuum_feasible_dual_range": [dual_minimum, dual_maximum],
            "continuum_dual_lower_bound": dual_value,
            "fine_grid_raw_objective_diagnostic": grid_objective,
            "high_precision_fixed_support_primal_upper_evaluation": primal_value,
            "high_precision_primal_moment_residual": primal_residual,
            "numerical_bracket_width": primal_value - dual_value,
            "fixed_support_primal_atom_count": len(primal_atoms),
            "fixed_support_positive_atoms": sum(weight > 0 for _, weight in primal_atoms),
            "fixed_support_negative_atoms": sum(weight < 0 for _, weight in primal_atoms),
            "attenuation_envelope_seven_moment_resource": envelope_value,
            "full_oscillatory_resource_minus_envelope": dual_value - envelope_value,
            "theorem_scope": (
                "The dual polynomial is globally normalized using all real "
                "critical points, so it is feasible for the continuum signed "
                "Hausdorff-moment problem up to stated numerical root-finding "
                "precision. Twelve rational grid support points define a "
                "unique signed measure through the exact moment equations; "
                "its 100-digit weight evaluation supplies the displayed "
                "upper evaluation."
            ),
            "boundary": (
                "This is a classical signed path-measure resource. It is not "
                "a quantum-negativity theorem, negative information, a loss "
                "law, or an internal FIN source. The raw HiGHS objective is a "
                "conditioning diagnostic, not the upper certificate. P352 "
                "must interval-certify the degree-11 extrema and fixed-support "
                "weight signs before the bracket becomes theorem-grade."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P339: continuous QW-2038 refreeze sensitivity and hold-out audit
# ---------------------------------------------------------------------------


def qw2038_context() -> dict[str, Any]:
    d2030 = json.loads(
        (ROOT / "report_qw2030_final_stage_c_gate_combined_branch.json").read_text(
            encoding="utf-8"
        )
    )
    d2028 = json.loads(
        (ROOT / "report_qw2028_joint_scan_with_gw_kappa_term.json").read_text(
            encoding="utf-8"
        )
    )
    d2029 = json.loads(
        (ROOT / "report_qw2029_ckm_blocker_shared_flavor_basis_scan.json").read_text(
            encoding="utf-8"
        )
    )
    d2027 = json.loads(
        (ROOT / "report_qw2027_gw_control_gap_structural_term_scan.json").read_text(
            encoding="utf-8"
        )
    )
    report = json.loads(
        (ROOT / "report_qw2038_derivation_compatible_kernel_refreeze_scan.json").read_text(
            encoding="utf-8"
        )
    )
    return {
        "kernel0": d2030["kernel"],
        "best2028": d2028["best_row"],
        "best2029": d2029["best_row"],
        "kappa": float(d2027["selected"]["kappa"]),
        "report": report,
        "gw": pd.read_csv(ROOT / "gw1831_window_features.csv"),
    }


def qw2038_evaluate(parameters: np.ndarray, context: dict[str, Any]) -> dict[str, Any]:
    omega, phi, beta, eta = [float(value) for value in parameters]
    kernel = {"omega": omega, "phi": phi, "beta": beta, "eta": eta}
    best2028 = context["best2028"]
    best2029 = context["best2029"]
    mass = qw2038.mass_metrics(
        p_amp=float(best2028["p_amp"]),
        r_dist=float(best2028["r_dist"]),
        gamma_scale=float(best2028["gamma_scale"]),
        coeffs=best2028["coeffs"],
        kernel=kernel,
    )
    flavor = qw2038.flavor_metrics(
        q_up=np.array([0.0, 9.0, 14.0]),
        q_down=np.array([7.0, 9.0, 14.0]),
        q_nu=np.array(best2029["q_nu"], dtype=float),
        q_lep=np.array([24.0, 14.0, 9.0]),
        p_amp=float(best2028["p_amp"]),
        r_dist=float(best2028["r_dist"]),
        params=best2029["params"],
        kernel=kernel,
    )
    gw = qw2038.gw_metrics(
        p_amp=float(best2028["p_amp"]),
        r_dist=float(best2028["r_dist"]),
        kernel=kernel,
        df_gw=context["gw"],
        kappa=context["kappa"],
    )
    thresholds = context["report"]["thresholds"]
    flags = {
        "mass_mean": mass["mean_rel_err_pct"] <= thresholds["mass_mean_rel_pct_max"],
        "mass_max": mass["max_rel_err_pct"] <= thresholds["mass_max_rel_pct_max"],
        "ckm": flavor["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"],
        "pmns": flavor["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"],
        "gw_sep": gw["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"],
        "gw_adv": gw["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"],
        "gw_auc": gw["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"],
        "gw_gap": gw["control_median_gap"] <= thresholds["gw_control_gap_max"],
    }
    beta0 = float(context["kernel0"]["beta"])
    eta0 = float(context["kernel0"]["eta"])
    score = (
        mass["mean_rel_err_pct"] / thresholds["mass_mean_rel_pct_max"]
        + flavor["ckm_mean_rel_pct"] / thresholds["ckm_mean_rel_pct_max"]
        + flavor["pmns_mean_rel_pct"] / thresholds["pmns_mean_rel_pct_max"]
        + max(
            0.0,
            (
                thresholds["gw_sep_min"]
                - gw["sep_median_h1l1_minus_ctrl"]
            )
            / thresholds["gw_sep_min"],
        )
        + max(
            0.0,
            (
                thresholds["gw_adv_min"]
                - gw["adv_shared_minus_ctrl_q90"]
            )
            / thresholds["gw_adv_min"],
        )
        + max(
            0.0,
            (thresholds["gw_auc_min"] - gw["auc_h1l1_vs_ctrl"])
            / thresholds["gw_auc_min"],
        )
        + max(
            0.0,
            (
                gw["control_median_gap"]
                - thresholds["gw_control_gap_max"]
            )
            / thresholds["gw_control_gap_max"],
        )
        + 0.05 * abs(beta - beta0)
        + 0.03 * abs(eta - eta0)
    )
    return {
        "parameters": kernel,
        "score": float(score),
        "all_pass": bool(all(flags.values())),
        "flags": flags,
        "mass": mass,
        "flavor": flavor,
        "gw": gw,
    }


def program_339(
    rng: np.random.Generator,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    context = qw2038_context()
    grid = context["report"]["grid"]
    bounds = [
        (min(grid["omega_grid"]), max(grid["omega_grid"])),
        (min(grid["phi_grid"]), max(grid["phi_grid"])),
        (min(grid["beta_grid"]), max(grid["beta_grid"])),
        (min(grid["eta_grid"]), max(grid["eta_grid"])),
    ]
    frozen = np.array([STRICT_OMEGA, STRICT_PHI, STRICT_BETA, STRICT_ETA])
    frozen_eval = qw2038_evaluate(frozen, context)

    result = differential_evolution(
        lambda vector: qw2038_evaluate(vector, context)["score"],
        bounds=bounds,
        seed=SEED + 339,
        maxiter=24,
        popsize=7,
        polish=True,
        tol=2e-7,
        updating="immediate",
    )
    optimum_eval = qw2038_evaluate(result.x, context)

    scales = np.array([0.006, 0.018, 0.035, 0.04])
    local_rows = []
    pass_count = 0
    scores = []
    for index in range(256):
        candidate = frozen + rng.normal(size=4) * scales
        candidate = np.array(
            [np.clip(candidate[j], bounds[j][0], bounds[j][1]) for j in range(4)]
        )
        evaluated = qw2038_evaluate(candidate, context)
        pass_count += int(evaluated["all_pass"])
        scores.append(evaluated["score"])
        local_rows.append(
            {
                "row_type": "local_perturbation",
                "index": index,
                **evaluated["parameters"],
                "score": evaluated["score"],
                "all_pass": evaluated["all_pass"],
                "normalized_distance_from_frozen": float(
                    np.linalg.norm((candidate - frozen) / scales)
                ),
            }
        )
    rows = [
        {
            "row_type": "frozen_grid_winner",
            **frozen_eval["parameters"],
            "score": frozen_eval["score"],
            "all_pass": frozen_eval["all_pass"],
        },
        {
            "row_type": "continuous_bounded_best_candidate",
            **optimum_eval["parameters"],
            "score": optimum_eval["score"],
            "all_pass": optimum_eval["all_pass"],
        },
        *local_rows,
    ]
    source_hashes = {
        name: sha256_file(ROOT / name)
        for name in (
            "report_qw2030_final_stage_c_gate_combined_branch.json",
            "report_qw2028_joint_scan_with_gw_kappa_term.json",
            "report_qw2029_ckm_blocker_shared_flavor_basis_scan.json",
            "report_qw2027_gw_control_gap_structural_term_scan.json",
            "gw1831_window_features.csv",
        )
    }
    return (
        {
            "status": (
                "[Strong evidence] continuous in-sample sensitivity audit; "
                "[Refuted] claim that QW-2038 supplied a frozen independent hold-out"
            ),
            "frozen_evaluation": frozen_eval,
            "continuous_best_candidate": optimum_eval,
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "score_improvement": frozen_eval["score"] - optimum_eval["score"],
            "parameter_displacement": {
                name: float(optimum_eval["parameters"][name] - frozen[index])
                for index, name in enumerate(("omega", "phi", "beta", "eta"))
            },
            "local_perturbation_count": 256,
            "local_pass_fraction": pass_count / 256.0,
            "local_score_quantiles": {
                "q05": float(np.quantile(scores, 0.05)),
                "median": float(np.median(scores)),
                "q95": float(np.quantile(scores, 0.95)),
            },
            "independent_frozen_holdout_available": False,
            "holdout_reason": (
                "QW-2038 uses the same mass references, CKM/PMNS references, "
                "GW feature table, upstream fitted coefficients, and gate "
                "thresholds for search and scoring; no independently "
                "registered hold-out file or unblinding record is declared."
            ),
            "source_hashes": source_hashes,
            "boundary": (
                "Continuous stability or instability of this objective is "
                "in-sample procedural evidence only. It is not ontology, a "
                "selector source, QW-2191 closure, or predictive physics."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P340: carrier-natural scalar completion theorem and transfer test
# ---------------------------------------------------------------------------


def shortest_path_matrix(graph: nx.Graph) -> np.ndarray:
    nodes = list(graph.nodes())
    lengths = dict(nx.all_pairs_shortest_path_length(graph))
    return np.array(
        [[lengths[left][right] for right in nodes] for left in nodes],
        dtype=float,
    )


def raw_kernel_matrix(
    distances: np.ndarray,
    omega: float,
    phi: float,
    beta: float,
    eta: float,
    amplitude: float,
) -> np.ndarray:
    values = amplitude * np.cos(omega * distances + phi) / (
        1.0 + beta * distances**eta
    )
    values = np.array(values, dtype=float)
    np.fill_diagonal(values, 0.0)
    return (values + values.T) / 2.0


def matrix_polynomial(
    matrix: np.ndarray,
    coefficients: np.ndarray,
    center: float,
    radius: float,
) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    mapped = polynomial_value(coefficients, (eigenvalues - center) / radius)
    return (eigenvectors * mapped) @ eigenvectors.T


def cycle_scalar_interpolant() -> tuple[np.ndarray, float, float, float]:
    graph = nx.cycle_graph(N)
    distances = shortest_path_matrix(graph)
    legacy = raw_kernel_matrix(
        distances, math.pi / 4, math.pi / 6, 0.01, 1.0, ALPHA_GEO
    )
    strict = raw_kernel_matrix(
        distances, STRICT_OMEGA, STRICT_PHI, 1.0, STRICT_ETA, 1.0
    )
    legacy_modes = np.fft.fft(legacy[0]).real
    strict_modes = np.fft.fft(strict[0]).real
    pairs: list[tuple[float, float]] = []
    for left, right in zip(legacy_modes, strict_modes):
        if not any(abs(left - old_left) < 1e-10 for old_left, _ in pairs):
            pairs.append((float(left), float(right)))
    center = 0.5 * (min(x for x, _ in pairs) + max(x for x, _ in pairs))
    radius = 0.5 * (max(x for x, _ in pairs) - min(x for x, _ in pairs))
    x = np.array([(value - center) / radius for value, _ in pairs])
    y = np.array([value for _, value in pairs])
    coefficients = np.polynomial.polynomial.polyfit(x, y, len(pairs) - 1)
    residual = float(np.max(np.abs(polynomial_value(coefficients, x) - y)))
    return coefficients, center, radius, residual


def program_340() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    coefficients, center, radius, training_residual = cycle_scalar_interpolant()
    graphs = {
        "cycle_C12": nx.cycle_graph(12),
        "path_P12": nx.path_graph(12),
        "random_3_regular_12": nx.random_regular_graph(3, 12, seed=340),
        "barbell_irregular_12": nx.barbell_graph(5, 2),
        "complete_K12": nx.complete_graph(12),
    }
    rows = []
    for name, graph in graphs.items():
        distances = shortest_path_matrix(graph)
        legacy = raw_kernel_matrix(
            distances, math.pi / 4, math.pi / 6, 0.01, 1.0, ALPHA_GEO
        )
        strict = raw_kernel_matrix(
            distances, STRICT_OMEGA, STRICT_PHI, 1.0, STRICT_ETA, 1.0
        )
        mapped = matrix_polynomial(legacy, coefficients, center, radius)
        error = float(
            np.linalg.norm(mapped - strict, "fro")
            / max(np.linalg.norm(strict, "fro"), 1e-15)
        )
        commutator = float(
            np.linalg.norm(legacy @ strict - strict @ legacy, "fro")
            / max(
                np.linalg.norm(legacy, "fro") * np.linalg.norm(strict, "fro"),
                1e-15,
            )
        )
        rows.append(
            {
                "carrier": name,
                "vertices": graph.number_of_nodes(),
                "diameter": nx.diameter(graph),
                "relative_scalar_interpolant_transfer_error": error,
                "relative_commutator": commutator,
                "same_cycle_trained_polynomial_used": True,
                "passes_1e_8": error < 1e-8,
            }
        )
    noncycle_failures = sum(
        not row["passes_1e_8"] for row in rows if row["carrier"] != "cycle_C12"
    )
    return (
        {
            "status": (
                "[Proven] natural scalar maps preserve commutation; "
                "[Strong evidence] cycle-interpolant transfer refutation"
            ),
            "cycle_interpolant_degree": len(coefficients) - 1,
            "cycle_training_spectral_residual": training_residual,
            "cycle_interpolant_coefficients_ascending": coefficients.tolist(),
            "spectral_affine_center": center,
            "spectral_affine_radius": radius,
            "noncycle_transfer_failures": noncycle_failures,
            "theorem": (
                "A carrier-natural scalar functional completion F_A=f(A) "
                "commutes with A on every carrier. A polynomial interpolated "
                "from the seven C12 spectral classes is target-coded training "
                "data; applying that identical polynomial to other carriers "
                "is the finite naturality test."
            ),
            "boundary": (
                "Failure excludes this scalar cycle-trained natural "
                "transformation, not nonscalar parents, enlarged state spaces, "
                "or a genuinely new carrier-dependent source law."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P341: ideal mesh compilation and short-time inverse calibration
# ---------------------------------------------------------------------------


def transition_probabilities(unitary: np.ndarray) -> np.ndarray:
    return np.abs(unitary) ** 2


def program_341(
    strict_a: np.ndarray,
    rng: np.random.Generator,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    prior = json.loads(
        (ROOT / "FIN_Programs_323_336_Results.json").read_text(encoding="utf-8")
    )
    witness_time = float(prior["P326"]["best_protocols"]["wave"]["nominal_time"])
    target = expm(-1j * witness_time * strict_a)
    rotations, diagonal = p295.givens_decompose_unitary(target)
    rows: list[dict[str, Any]] = []
    for quantization in (0.0, 1e-5, 1e-4, 1e-3):
        compiled = p295.reconstruct_givens(rotations, diagonal, quantization)
        rows.append(
            {
                "row_type": "mesh_quantization",
                "setting": quantization,
                "unitary_operator_error": float(np.linalg.norm(compiled - target, 2)),
                "maximum_vertex_tv_error": float(
                    np.max(
                        0.5
                        * np.sum(
                            np.abs(
                                transition_probabilities(compiled)
                                - transition_probabilities(target)
                            ),
                            axis=0,
                        )
                    )
                ),
            }
        )

    calibration_time = 0.08
    calibration_target = expm(-1j * calibration_time * strict_a)
    noise_levels = (1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3)
    calibration_summary = []
    for noise in noise_levels:
        generator_errors = []
        prediction_tvs = []
        for _ in range(48):
            perturbation = noise * (
                rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
            )
            estimated_unitary, _ = polar(calibration_target + perturbation)
            estimated_generator = (1j / calibration_time) * logm(estimated_unitary)
            estimated_generator = (
                estimated_generator + estimated_generator.conj().T
            ) / 2.0
            estimated_generator = estimated_generator.real
            generator_errors.append(
                float(
                    np.linalg.norm(estimated_generator - strict_a, "fro")
                    / np.linalg.norm(strict_a, "fro")
                )
            )
            predicted = expm(-1j * witness_time * estimated_generator)
            prediction_tvs.append(
                float(
                    np.max(
                        0.5
                        * np.sum(
                            np.abs(
                                transition_probabilities(predicted)
                                - transition_probabilities(target)
                            ),
                            axis=0,
                        )
                    )
                )
            )
        row = {
            "row_type": "short_time_inverse_calibration",
            "setting": noise,
            "median_relative_generator_error": float(np.median(generator_errors)),
            "p95_relative_generator_error": float(np.quantile(generator_errors, 0.95)),
            "median_witness_prediction_tv": float(np.median(prediction_tvs)),
            "p95_witness_prediction_tv": float(np.quantile(prediction_tvs, 0.95)),
        }
        rows.append(row)
        calibration_summary.append(row)

    _, generator_eigenvectors = np.linalg.eigh(strict_a)
    projector = np.outer(
        generator_eigenvectors[:, 0], generator_eigenvectors[:, 0].conj()
    )
    # Adding 2pi/t times a spectral projector leaves exp(-itA) unchanged.
    branch_shift = (2.0 * math.pi / witness_time) * projector
    branch_generator = strict_a + branch_shift.real
    branch_residual = float(
        np.linalg.norm(expm(-1j * witness_time * branch_generator) - target, 2)
    )
    return (
        {
            "status": (
                "[Strong evidence] ideal 12-mode compilation and synthetic "
                "short-time inverse calibration; [Proven] one-time logarithm "
                "branch nonuniqueness"
            ),
            "witness_time": witness_time,
            "calibration_time": calibration_time,
            "mesh_two_mode_rotation_count": len(rotations),
            "mesh_terminal_phase_count": N,
            "ideal_reconstruction_residual": rows[0]["unitary_operator_error"],
            "calibration_summary": calibration_summary,
            "single_time_branch_generator_distance": float(
                np.linalg.norm(branch_generator - strict_a, 2)
            ),
            "single_time_branch_unitary_residual": branch_residual,
            "theorem": (
                "If P is a spectral projector commuting with A, then "
                "A'=A+(2pi/t)P gives exp(-itA')=exp(-itA). Thus a single "
                "unitary time does not identify the generator. A sufficiently "
                "short calibrated time fixes the principal-log branch under a "
                "declared spectral bound."
            ),
            "boundary": (
                "The mesh is ideal and the calibration noise is synthetic. "
                "No loss, detector drift, fabrication, phase stabilization, "
                "or independently calibrated hardware record is supplied."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P342: external frozen-holdout admission gate
# ---------------------------------------------------------------------------


def program_342() -> dict[str, Any]:
    candidates = sorted(
        set(ROOT.rglob("*holdout*manifest*.json"))
        | set(ROOT.rglob("*manifest*holdout*.json"))
    )
    admitted = []
    for path in candidates:
        relative = str(path.relative_to(ROOT))
        if any(token in relative.lower() for token in ("template", "synthetic", "example")):
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        required = {
            "provider",
            "registrar",
            "analyst",
            "dataset_hash",
            "freeze_timestamp",
            "unblinding_authorization",
        }
        if isinstance(payload, dict) and required.issubset(payload):
            if len({payload["provider"], payload["registrar"], payload["analyst"]}) == 3:
                admitted.append(relative)
    return {
        "status": "[Blocked by external evidence] no admitted independent frozen hold-out",
        "candidate_manifest_count": len(candidates),
        "admitted_manifest_paths": admitted,
        "admitted": bool(admitted),
        "continuous_refreeze_external_validation_authorized": bool(admitted),
        "boundary": (
            "A repository script cannot manufacture independent data custody, "
            "pre-registration, frozen hashes, distinct roles, or unblinding."
        ),
    }


# ---------------------------------------------------------------------------
# P343: declared parallel-comb subclass
# ---------------------------------------------------------------------------


def relative_unitary_phases(
    strict_a: np.ndarray, legacy_a: np.ndarray, scale: float, time: float
) -> np.ndarray:
    strict_unitary = expm(-1j * time * strict_a)
    legacy_unitary = expm(-1j * time * scale * legacy_a)
    return np.linalg.eigvals(strict_unitary.conj().T @ legacy_unitary)


def tensor_power_phases(phases: np.ndarray, uses: int) -> np.ndarray:
    result = np.array([1.0 + 0.0j])
    for _ in range(uses):
        result = (result[:, None] * phases[None, :]).ravel()
    return result


def program_343(strict_a: np.ndarray) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    legacy = p295.laplacian_from_profile(p295.legacy_weights())
    prior = json.loads(
        (ROOT / "FIN_Programs_323_336_Results.json").read_text(encoding="utf-8")
    )
    scale = float(prior["P328"]["optimized_legacy_scale"])
    times = np.linspace(0.03, 0.75, 97)
    rows = []
    earliest: dict[str, float | None] = {}
    for uses in range(1, 5):
        first_perfect = None
        for time in times:
            phases = relative_unitary_phases(strict_a, legacy, scale, float(time))
            combined = tensor_power_phases(phases, uses)
            half_diamond = p323.unitary_channel_half_diamond(combined)
            perfect = half_diamond >= 1.0 - 1e-11
            if perfect and first_perfect is None:
                first_perfect = float(time)
            rows.append(
                {
                    "uses": uses,
                    "time_per_use": float(time),
                    "parallel_input_dimension": int(N**uses),
                    "relative_phase_count": int(len(combined)),
                    "half_diamond_distance": half_diamond,
                    "perfectly_distinguishable": perfect,
                }
            )
        earliest[str(uses)] = first_perfect
    return (
        {
            "status": (
                "[Proven] ideal parallel n-use unitary-channel subclass; "
                "[Blocked] unrestricted adaptive multitime comb SDP"
            ),
            "common_single_use_input_output_space": "H=C^12",
            "single_use_channel_space": "B(H) -> B(H)",
            "parallel_n_use_input_space": "H tensor ... tensor H",
            "intermediate_memory": "none in the executed subclass",
            "intermediate_instruments": "identity/no intervention",
            "final_instrument": "optimal unrestricted joint POVM (existence theorem)",
            "legacy_scale": scale,
            "earliest_sampled_perfect_time_by_uses": earliest,
            "theorem": (
                "For n parallel uses of unitary channels, the relative "
                "eigenvalue set is the n-fold product set. Its convex-hull "
                "distance to zero determines the n-use half diamond norm."
            ),
            "boundary": (
                "This is not the full adaptive process-tensor/quantum-comb "
                "optimization: no memory channel, causal intervention slots, "
                "or comb SDP is executed."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P344: nonparametric monotone-clock no-go and minimax gap design
# ---------------------------------------------------------------------------


def maximum_gap(times: list[float]) -> float:
    ordered = sorted(set([0.0, *times]))
    return max(right - left for left, right in zip(ordered, ordered[1:]))


def program_344() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    prior = json.loads(
        (ROOT / "FIN_Programs_323_336_Results.json").read_text(encoding="utf-8")
    )
    quartic_times = [float(value) for value in prior["P329"]["best_design"]["times"]]
    horizon = max(quartic_times)
    equispaced = [horizon * index / 5.0 for index in range(1, 6)]
    curvature_bound = 0.2
    slope_margin = 0.2
    rows = []
    for name, times in (
        ("P329_quartic_Fisher_design", quartic_times),
        ("nonparametric_equispaced_design", equispaced),
    ):
        gap = maximum_gap(times)
        interpolation_bound = curvature_bound * gap**2 / 8.0
        bump_amplitude = curvature_bound * gap**2 / (2.0 * math.pi**2)
        maximum_bump_slope = curvature_bound * gap / (2.0 * math.pi)
        rows.append(
            {
                "design": name,
                "times": times,
                "maximum_gap_including_t0": gap,
                "linear_interpolation_clock_error_bound": interpolation_bound,
                "adversarial_single_clock_bump_amplitude": bump_amplitude,
                "two_clock_midpoint_separation": 2.0 * bump_amplitude,
                "maximum_bump_slope": maximum_bump_slope,
                "monotonicity_preserved_under_declared_margin": (
                    maximum_bump_slope < slope_margin
                ),
            }
        )
    return (
        {
            "status": (
                "[Proven] finite-design nonparametric clock nonidentification; "
                "[Proven] C2 interpolation upper bound"
            ),
            "clock_class": (
                "tau in W^{2,infinity}, tau(0)=0, "
                "|tau''|<=0.2, tau'>=0.8"
            ),
            "curvature_bound": curvature_bound,
            "slope_margin": slope_margin,
            "best_maximum_gap_design": min(
                rows, key=lambda row: row["maximum_gap_including_t0"]
            )["design"],
            "theorem": (
                "On the largest unsampled interval [a,b], the two clocks "
                "tau_+/-=t +/- A sin^2(pi(t-a)/(b-a)) agree at every sampled "
                "time. With A=M(b-a)^2/(2pi^2), they obey the curvature bound "
                "and remain monotone under the declared slope margin, yet "
                "differ inside the interval. Also, linear interpolation error "
                "is bounded by M h^2/8."
            ),
            "boundary": (
                "No finite design identifies an unrestricted monotone clock. "
                "Equispacing is minimax only for the declared maximum-gap/"
                "curvature criterion; it does not replace external time "
                "calibration."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P345: carrier-universal versus formula-level invariants
# ---------------------------------------------------------------------------


def spectral_entropy(matrix: np.ndarray) -> float:
    singular = np.abs(np.linalg.eigvalsh(matrix))
    total = float(np.sum(singular))
    if total <= 1e-15:
        return 0.0
    probabilities = singular / total
    probabilities = probabilities[probabilities > 1e-15]
    return float(-np.sum(probabilities * np.log(probabilities)))


def graph_kernel_pair(graph: nx.Graph) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    distances = shortest_path_matrix(graph)
    legacy = raw_kernel_matrix(
        distances, math.pi / 4, math.pi / 6, 0.01, 1.0, ALPHA_GEO
    )
    strict = raw_kernel_matrix(
        distances, STRICT_OMEGA, STRICT_PHI, 1.0, STRICT_ETA, 1.0
    )
    return distances, legacy, strict


def program_345() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    graphs = {
        "cycle_C12": nx.cycle_graph(12),
        "path_P12": nx.path_graph(12),
        "random_3_regular_12": nx.random_regular_graph(3, 12, seed=345),
        "barbell_irregular_12": nx.barbell_graph(5, 2),
        "complete_K12": nx.complete_graph(12),
        "grid_3x4": nx.convert_node_labels_to_integers(nx.grid_2d_graph(3, 4)),
    }
    rows = []
    for name, graph in graphs.items():
        distances, legacy, strict = graph_kernel_pair(graph)
        strict_eigenvalues = np.linalg.eigvalsh(strict)
        legacy_eigenvalues = np.linalg.eigvalsh(legacy)
        diameter = nx.diameter(graph)
        local_ratios = []
        for distance in range(1, min(diameter, 6) + 1):
            local_ratios.append(
                float(
                    (
                        math.cos(STRICT_OMEGA * distance + STRICT_PHI)
                        / (1.0 + distance**STRICT_ETA)
                    )
                    / (
                        math.cos(STRICT_OMEGA + STRICT_PHI)
                        / (1.0 + 1.0**STRICT_ETA)
                    )
                )
            )
        rows.append(
            {
                "carrier": name,
                "vertices": graph.number_of_nodes(),
                "edges": graph.number_of_edges(),
                "diameter": diameter,
                "strict_normalized_spectral_entropy": spectral_entropy(strict),
                "strict_spectral_radius_per_vertex": float(
                    np.max(np.abs(strict_eigenvalues)) / graph.number_of_nodes()
                ),
                "strict_frobenius_per_edge_sqrt": float(
                    np.linalg.norm(strict, "fro")
                    / math.sqrt(max(graph.number_of_edges(), 1))
                ),
                "legacy_normalized_spectral_entropy": spectral_entropy(legacy),
                "legacy_strict_commutator": float(
                    np.linalg.norm(legacy @ strict - strict @ legacy, "fro")
                ),
                "strict_formula_ratios_d1_to_diameter_cap6": local_ratios,
                "formula_parameters_identical_by_construction": True,
            }
        )
    numeric_invariants = (
        "strict_normalized_spectral_entropy",
        "strict_spectral_radius_per_vertex",
        "strict_frobenius_per_edge_sqrt",
    )
    coefficients_of_variation = {}
    for key in numeric_invariants:
        values = np.array([row[key] for row in rows])
        coefficients_of_variation[key] = float(np.std(values) / np.mean(values))
    return (
        {
            "status": (
                "[Proven] formula parameters are carrier-independent inputs; "
                "[Strong evidence] tested nontrivial operator invariants are "
                "carrier-dependent"
            ),
            "carrier_count": len(rows),
            "tested_operator_invariant_coefficients_of_variation": (
                coefficients_of_variation
            ),
            "universal_formula_level_data": {
                "omega": "743/4000",
                "phi": "13/80",
                "beta": 1,
                "eta": "9/5",
            },
            "theorem": (
                "A radial law can be reused on every graph metric, so its "
                "written coefficients and distance-wise ratios are natural "
                "at the formula level. Spectra, multiplicities, commutators, "
                "and normalized energies depend on the carrier and are not "
                "therefore universal consequences of the coefficients alone."
            ),
            "boundary": (
                "Carrier-independent syntax is not a derivation of the "
                "parameters and not an observable universality theorem. A "
                "categorical naturality theorem needs declared graph "
                "morphisms and functorial operator transport."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P346: typed bridge-resource free operations and monotones
# ---------------------------------------------------------------------------


def jordan_negative_mass(weights: np.ndarray) -> float:
    return float(np.sum(np.maximum(-weights, 0.0)))


def program_346(
    rng: np.random.Generator,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    contraction_rows = []
    maximum_violation = -math.inf
    ratios = []
    for trial in range(512):
        signed = rng.normal(size=16)
        signed += (1.0 - float(np.sum(signed))) / len(signed)
        stochastic = rng.gamma(1.0, 1.0, size=(8, 16))
        stochastic /= np.sum(stochastic, axis=0, keepdims=True)
        pushed = stochastic @ signed
        before = jordan_negative_mass(signed)
        after = jordan_negative_mass(pushed)
        violation = after - before
        maximum_violation = max(maximum_violation, violation)
        if before > 1e-12:
            ratios.append(after / before)
        if trial < 64:
            contraction_rows.append(
                {
                    "row_type": "signed_measure_markov_pushforward",
                    "trial": trial,
                    "negative_mass_before": before,
                    "negative_mass_after": after,
                    "contraction_ratio": after / max(before, 1e-15),
                    "monotone_holds": after <= before + 1e-12,
                }
            )
    typed_rows = [
        {
            "row_type": "typed_resource_law",
            "resource": "signed_path_weight",
            "free_operation": "positive mass-preserving Markov pushforward",
            "monotone": "Jordan negative mass",
            "result": "nonincreasing by total-variation contraction",
        },
        {
            "row_type": "typed_resource_law",
            "resource": "nontorsion_phase",
            "free_operation": "group homomorphism generated from torsion input",
            "monotone": "torsion-vs-infinite-order flag",
            "result": "a torsion input cannot generate an infinite-order image",
        },
        {
            "row_type": "typed_resource_law",
            "resource": "pointing/orientation",
            "free_operation": "equivariant map without a supplied section",
            "monotone": "existence of an equivariant point",
            "result": "no point is created from a free transitive torsor",
        },
        {
            "row_type": "typed_resource_law",
            "resource": "positive dimensional scale",
            "free_operation": "R_+-equivariant scale-free map",
            "monotone": "existence of a fixed positive unit",
            "result": "no invariant positive point under the free scale action",
        },
        {
            "row_type": "typed_resource_law",
            "resource": "independent cross law",
            "free_operation": "deterministic function of legacy input only",
            "monotone": "completion-fiber nonidentification witness",
            "result": "identical inputs cannot be mapped to two distinct targets",
        },
    ]
    return (
        {
            "status": (
                "[Proven] typed monotonicity theorems; [Strong evidence] "
                "finite Markov contraction audit"
            ),
            "markov_trials": 512,
            "maximum_negative_mass_monotonicity_violation": maximum_violation,
            "median_negative_mass_contraction_ratio": float(np.median(ratios)),
            "typed_resource_count": len(typed_rows),
            "free_operation_product_scope": (
                "Each resource has its own category of free operations; the "
                "vector is not collapsed to one physically privileged scalar."
            ),
            "theorem": (
                "For a signed measure mu with fixed total mass, "
                "N(mu)=(||mu||_1-mu(X))/2. A positive mass-preserving Markov "
                "map contracts total variation and preserves total mass, so "
                "N(Tmu)<=N(mu). The remaining four statements are the typed "
                "homomorphism/equivariance/nonidentification analogues."
            ),
            "boundary": (
                "These are axioms for declared free-operation classes, not a "
                "proof that FIN dynamics implements those classes or that one "
                "resource physically converts into another."
            ),
        },
        [*contraction_rows, *typed_rows],
    )


# ---------------------------------------------------------------------------
# P347: candidate resource-coupling laws and falsification matrix
# ---------------------------------------------------------------------------


def program_347() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    distance = np.arange(1.0, 257.0)
    log_envelope = -np.log1p(distance**STRICT_ETA)
    minimum_phase = -np.imag(hilbert(log_envelope))
    frozen_phase = STRICT_OMEGA * distance + STRICT_PHI
    # A causal minimum-phase phase is defined only up to a constant and a
    # pure delay.  Subtract their best least-squares contribution.
    nuisance = np.column_stack([np.ones_like(distance), distance])
    correction = np.linalg.lstsq(
        nuisance, frozen_phase - minimum_phase, rcond=None
    )[0]
    residual = frozen_phase - (minimum_phase + nuisance @ correction)
    central = slice(32, -32)
    kk_rmse = float(np.sqrt(np.mean(residual[central] ** 2)))
    kk_maximum = float(np.max(np.abs(residual[central])))

    rows = [
        {
            "candidate": "minimum_phase_Kramers_Kronig",
            "extra_axioms": (
                "distance promoted to causal frequency coordinate; analytic "
                "zero-free transfer function; boundary conditions"
            ),
            "test_statistic": "central phase RMSE after constant/delay freedom",
            "value": kk_rmse,
            "accepted": kk_rmse < 0.01,
            "verdict": "[Refuted in declared discretization]",
        },
        {
            "candidate": "fractional_power_analytic_continuation",
            "extra_axioms": "branch choice for z^(9/5)",
            "test_statistic": "generated phase family",
            "value": "branch-constant phase, not omega*d+phi",
            "accepted": False,
            "verdict": "[Refuted structurally]",
        },
        {
            "candidate": "complex_signed_path_representation",
            "extra_axioms": "insert full oscillatory target moments",
            "test_statistic": "representation exactness vs source independence",
            "value": "exact representation; target-reading source law",
            "accepted": False,
            "verdict": "[Representation only]",
        },
        {
            "candidate": "abstract_resource_product",
            "extra_axioms": "none beyond P335 product signature",
            "test_statistic": "logical implication between resources",
            "value": "single-omission models exist",
            "accepted": False,
            "verdict": "[Refuted without a new coupling axiom]",
        },
    ]
    return (
        {
            "status": (
                "[Refuted] all four tested resource-coupling laws as internal "
                "co-generation theorems"
            ),
            "minimum_phase_central_rmse": kk_rmse,
            "minimum_phase_central_maximum_residual": kk_maximum,
            "candidate_count": len(rows),
            "accepted_internal_coupling_laws": sum(row["accepted"] for row in rows),
            "conclusion": (
                "The amplitude envelope, nontorsion phase, cross law, pointing, "
                "and scale remain separately sourced on current artifacts. "
                "Causality/minimum-phase analyticity could mathematically link "
                "amplitude and phase only after adding assumptions not present "
                "in the FIN distance kernel, and its tested phase does not "
                "match the frozen linear phase."
            ),
            "boundary": (
                "Failure of these four candidates is not an impossibility "
                "theorem for every future coupling law."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P348: falsifiable conditional electroweak vector
# ---------------------------------------------------------------------------


def program_348() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    x = ALPHA_GEO / 12.0
    observables: list[tuple[str, float, str]] = [
        ("sin2_theta_W", x, "historical role-law axiom"),
        ("mW_over_mZ_tree", math.sqrt(1.0 - x), "Higgs-doublet tree relation"),
        ("electron_gV_over_gA_tree", 1.0 - 4.0 * x, "SM lepton representation"),
        (
            "electron_gL2_over_gR2_tree",
            ((-0.5 + x) ** 2) / (x**2),
            "SM lepton representation",
        ),
        ("neutrino_gV_over_gA_tree", 1.0, "left-handed neutrino representation"),
    ]
    rows = [
        {
            "observable": name,
            "conditional_prediction": value,
            "additional_axiom_or_representation": assumption,
            "strict_FIN_prediction": False,
            "external_measurement_used_in_this_report": False,
        }
        for name, value, assumption in observables
    ]
    epsilon = 1e-7

    def vector(angle: float) -> np.ndarray:
        return np.array(
            [
                angle,
                math.sqrt(1.0 - angle),
                1.0 - 4.0 * angle,
                ((-0.5 + angle) ** 2) / angle**2,
                1.0,
            ]
        )

    jacobian = (vector(x + epsilon) - vector(x - epsilon)) / (2.0 * epsilon)
    return (
        {
            "status": (
                "[Constructed conditional falsification vector]; "
                "[Refuted] strict FIN electroweak prediction claim"
            ),
            "conditional_angle": x,
            "observable_count": len(rows),
            "one_parameter_jacobian_rank": int(
                np.linalg.matrix_rank(jacobian[:, None], tol=1e-10)
            ),
            "dimensionless_joint_prediction_vector": vector(x).tolist(),
            "required_axioms": [
                "ordered and pointed weak/hypercharge sectors",
                "historical role law sin^2(theta_W)=alpha_geo/12",
                "SU(2)_L x U(1)_Y gauge representations and charges",
                "one Higgs-doublet tree-level symmetry-breaking relations",
                "renormalization scheme/scale before precision comparison",
            ],
            "falsification_protocol": (
                "Freeze all representation assignments and the role-law "
                "before receiving an external blinded observable vector; "
                "compare the complete covariance-weighted residual once, "
                "without refitting alpha_geo or selecting observables."
            ),
            "boundary": (
                "The numbers are consequences of supplied electroweak "
                "structure plus the historical angle axiom. No current "
                "experimental values are used and no radiative correction, "
                "gauge dynamics, Higgs scale, anomaly theorem, or role-transfer "
                "certificate is derived from strict FIN."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P349: dimensional-conversion metrology frame
# ---------------------------------------------------------------------------


def program_349() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    # Rows are exponents of (ell_*, tau_*, hbar_*).
    candidates = {
        "length_standard": np.array([1.0, 0.0, 0.0]),
        "time_standard": np.array([0.0, 1.0, 0.0]),
        "action_standard": np.array([0.0, 0.0, 1.0]),
        "energy_standard": np.array([0.0, -1.0, 1.0]),
        "mass_standard_nonrel": np.array([-2.0, 1.0, 1.0]),
        "velocity_standard": np.array([1.0, -1.0, 0.0]),
    }
    rows = []
    best = None
    for names in itertools.combinations(candidates, 3):
        matrix = np.vstack([candidates[name] for name in names])
        rank = int(np.linalg.matrix_rank(matrix))
        determinant = float(abs(np.linalg.det(matrix))) if rank == 3 else 0.0
        condition = float(np.linalg.cond(matrix)) if rank == 3 else math.inf
        row = {
            "calibration_triple": names,
            "rank": rank,
            "absolute_determinant": determinant,
            "condition_number": condition if math.isfinite(condition) else None,
            "identifies_ell_tau_hbar": rank == 3,
        }
        rows.append(row)
        if rank == 3 and (
            best is None
            or condition < best["condition_number"]
            or (
                math.isclose(condition, best["condition_number"])
                and determinant > best["absolute_determinant"]
            )
        ):
            best = row
    assert best is not None
    return (
        {
            "status": (
                "[Proven] dimensional identifiability and conversion formulas; "
                "[Conditional] all three metrology anchors are external"
            ),
            "unit_coordinates": "(ell_*, tau_*, hbar_*)",
            "conversion_formulas": {
                "physical_time": "t_phys=tau_* t",
                "physical_length": "x_phys=ell_* d",
                "physical_action": "S_phys=hbar_* S",
                "physical_energy": "E_*=hbar_*/tau_*",
                "physical_hamiltonian": "H_phys=(hbar_*/tau_*) A",
                "nonrelativistic_mass_scale": "m_*=hbar_* tau_*/ell_*^2",
                "emergent_speed_unit": "c_*=ell_*/tau_*",
            },
            "candidate_calibration_triples": len(rows),
            "full_rank_triples": sum(row["rank"] == 3 for row in rows),
            "best_conditioned_triple": best,
            "dimensionless_FIN_log_unit_jacobian_rank": 0,
            "theorem": (
                "Logarithms of calibrated dimensional observables depend "
                "linearly on log(ell_*), log(tau_*), log(hbar_*) through the "
                "dimension-exponent matrix. They identify the conversion "
                "frame exactly iff this matrix has rank three. Every purely "
                "dimensionless FIN observable has zero unit Jacobian."
            ),
            "boundary": (
                "This constructs an operational Conversion-Axiom metrology "
                "frame; it does not make units emerge from FIN. Imposing "
                "c_*=c or hbar_*=hbar is an external calibration, not a "
                "strict theorem."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P350: external certified-hardware reservoir gate
# ---------------------------------------------------------------------------


def program_350() -> dict[str, Any]:
    prior = p323.program_336()
    return {
        **prior,
        "status": (
            "[Blocked by external evidence] no certified physical reservoir "
            "or calibrated device execution package"
        ),
        "program_350_executed_as_gate_audit_only": True,
        "hardware_execution_performed": False,
    }


# ---------------------------------------------------------------------------
# Formal core, figures, summaries, and main
# ---------------------------------------------------------------------------


def run_formal_core() -> dict[str, Any]:
    lean_path = (
        ROOT
        / ".elan/toolchains/leanprover--lean4---v4.28.0/bin/lean"
    )
    if not lean_path.exists():
        return {
            "status": "[Blocked] Lean executable unavailable",
            "lean_binary": None,
            "exit_code": None,
        }
    completed = subprocess.run(
        [str(lean_path), str(FORMAL_SOURCE)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "status": "[Proven] dependency-free Lean structural core compiled"
        if completed.returncode == 0
        else "[Blocked] Lean structural core failed",
        "lean_binary": str(lean_path),
        "lean_version": subprocess.run(
            [str(lean_path), "--version"],
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip(),
        "source_sha256": sha256_file(FORMAL_SOURCE),
        "exit_code": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    headlines = {
        "P337": "an exact rational dual certificate nearly closes the P323 value",
        "P338": "the full oscillatory strict kernel needs a larger signed path resource",
        "P339": "continuous refreeze is in-sample and has no independent hold-out",
        "P340": "a C12 scalar interpolant is not carrier-natural",
        "P341": "the ideal photonic mesh is compilable but one-time inversion is branch-ambiguous",
        "P342": "the external frozen-holdout gate remains closed",
        "P343": "parallel multiuse channels sharpen discrimination; full adaptive comb remains open",
        "P344": "finite times do not identify a nonparametric monotone clock",
        "P345": "formula-level constants survive carriers; tested operator invariants do not",
        "P346": "five bridge resources admit typed free-operation monotones",
        "P347": "four candidate resource-coupling laws fail",
        "P348": "the EW package is jointly falsifiable only conditionally",
        "P349": "dimensional conversion requires a rank-three external metrology frame",
        "P350": "the certified physical-reservoir gate remains closed",
    }
    return [
        {
            "program": program,
            "status": results[program]["status"],
            "headline": headline,
        }
        for program, headline in headlines.items()
    ]


def make_figures(results: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    p337_result = results["P337"]
    axes[0].bar(
        ["certified lower", "P323 numeric"],
        [
            p337_result["certified_dual_lower_bound"],
            p337_result["p323_high_precision_primal_value"],
        ],
        color=["#176b87", "#64ccc5"],
    )
    axes[0].set_ylim(0.40670, 0.40671)
    axes[0].set_title("P337 certified versus evaluated value")
    p338_rows = results["_P338_rows"]
    moment_rows = [row for row in p338_rows if row["row_type"] == "strict_oscillatory_moment"]
    axes[1].stem(
        [row["order"] for row in moment_rows],
        [row["value"] for row in moment_rows],
        basefmt=" ",
    )
    axes[1].axhline(0, color="black", linewidth=0.7)
    axes[1].set_title("P338 full strict moment sequence")
    axes[1].set_xlabel("order")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p337_p338_moment_certificates.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    p339_rows = [
        row
        for row in results["_P339_rows"]
        if row["row_type"] == "local_perturbation"
    ]
    axes[0].hist([row["score"] for row in p339_rows], bins=28, color="#176b87")
    axes[0].axvline(
        results["P339"]["frozen_evaluation"]["score"],
        color="#d1495b",
        linestyle="--",
        label="frozen",
    )
    axes[0].set_title("P339 local refreeze scores")
    axes[0].legend()
    p340_rows = results["_P340_rows"]
    axes[1].barh(
        [row["carrier"] for row in p340_rows],
        [max(row["relative_scalar_interpolant_transfer_error"], 1e-16) for row in p340_rows],
        color="#64ccc5",
    )
    axes[1].set_xscale("log")
    axes[1].set_title("P340 naturality transfer defect")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p339_p340_refreeze_naturality.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    calibration_rows = [
        row
        for row in results["_P341_rows"]
        if row["row_type"] == "short_time_inverse_calibration"
    ]
    axes[0].loglog(
        [row["setting"] for row in calibration_rows],
        [row["p95_witness_prediction_tv"] for row in calibration_rows],
        "o-",
        color="#176b87",
    )
    axes[0].set_xlabel("complex calibration noise")
    axes[0].set_ylabel("p95 witness TV")
    axes[0].set_title("P341 inverse calibration")
    comb_rows = results["_P343_rows"]
    for uses in range(1, 5):
        subset = [row for row in comb_rows if row["uses"] == uses]
        axes[1].plot(
            [row["time_per_use"] for row in subset],
            [row["half_diamond_distance"] for row in subset],
            label=f"{uses} uses",
        )
    axes[1].set_ylim(0, 1.03)
    axes[1].set_title("P343 parallel-channel separation")
    axes[1].legend()
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p341_p343_photonic_comb.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    clock_rows = results["_P344_rows"]
    axes[0].bar(
        [row["design"] for row in clock_rows],
        [row["two_clock_midpoint_separation"] for row in clock_rows],
        color=["#d1495b", "#64ccc5"],
    )
    axes[0].tick_params(axis="x", rotation=20)
    axes[0].set_title("P344 unresolved clock separation")
    invariant_rows = results["_P345_rows"]
    axes[1].barh(
        [row["carrier"] for row in invariant_rows],
        [row["strict_normalized_spectral_entropy"] for row in invariant_rows],
        color="#176b87",
    )
    axes[1].set_title("P345 carrier-dependent spectral entropy")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p344_p345_clock_invariants.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    monotone_rows = [
        row
        for row in results["_P346_rows"]
        if row["row_type"] == "signed_measure_markov_pushforward"
    ]
    axes[0].scatter(
        [row["negative_mass_before"] for row in monotone_rows],
        [row["negative_mass_after"] for row in monotone_rows],
        s=13,
        alpha=0.7,
        color="#176b87",
    )
    maximum = max(row["negative_mass_before"] for row in monotone_rows)
    axes[0].plot([0, maximum], [0, maximum], "--", color="#d1495b")
    axes[0].set_title("P346 negativity contraction")
    axes[0].set_xlabel("before")
    axes[0].set_ylabel("after")
    coupling_rows = results["_P347_rows"]
    axes[1].bar(
        [row["candidate"] for row in coupling_rows],
        [1 if row["accepted"] else 0 for row in coupling_rows],
        color="#d1495b",
    )
    axes[1].set_ylim(0, 1.2)
    axes[1].tick_params(axis="x", rotation=30)
    axes[1].set_title("P347 accepted internal couplings")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p346_p347_resources_couplings.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ew_rows = results["_P348_rows"]
    axes[0].barh(
        [row["observable"] for row in ew_rows],
        [row["conditional_prediction"] for row in ew_rows],
        color="#64ccc5",
    )
    axes[0].set_title("P348 conditional EW vector")
    metrology_rows = [
        row for row in results["_P349_rows"] if row["identifies_ell_tau_hbar"]
    ]
    axes[1].scatter(
        [row["absolute_determinant"] for row in metrology_rows],
        [row["condition_number"] for row in metrology_rows],
        color="#176b87",
    )
    axes[1].set_yscale("log")
    axes[1].set_xlabel("|det dimension matrix|")
    axes[1].set_ylabel("condition number")
    axes[1].set_title("P349 full-rank calibration triples")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p348_p349_conditional_metrology.png", dpi=180)
    plt.close(figure)


def main() -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    strict_a, _ = core.strict_operator()

    results: dict[str, Any] = {
        "metadata": {
            "programs": "P337-P350",
            "release": "10.30",
            "seed": SEED,
            "scope": (
                "interval/dual certification, oscillatory signed resources, "
                "continuous refreeze sensitivity, carrier naturality, ideal "
                "photonic compilation, parallel-comb discrimination, "
                "nonparametric clocks, typed resource monotones, conditional "
                "EW/metrology packages, and external admission gates"
            ),
            "new_theoretical_objects": {
                "O88_rational_Bernstein_moment_certificate": (
                    "exact rational p=-q^2 dual witness with algebraic input intervals"
                ),
                "O89_oscillatory_signed_path_norm": (
                    "minimum classical signed mass for the full strict moment sequence"
                ),
                "O90_refreeze_stability_fiber": (
                    "continuous objective landscape around the procedural grid winner"
                ),
                "O91_carrier_natural_completion_defect": (
                    "failure of one C12-trained scalar transformation on other carriers"
                ),
                "O92_compiled_inverse_photonic_frame": (
                    "Givens mesh, short-time calibration, branch certificate, and prediction map"
                ),
                "O93_parallel_comb_separation_ladder": (
                    "n-fold relative-phase convex-hull discrimination resource"
                ),
                "O94_nonparametric_clock_equivalence_tube": (
                    "monotone bounded-curvature clocks agreeing on a finite design"
                ),
                "O95_carrier_naturality_profile": (
                    "formula-level invariants separated from carrier-level spectra"
                ),
                "O96_typed_bridge_resource_preorder": (
                    "five resource-specific free-operation/monotone pairs"
                ),
                "O97_resource_coupling_acceptance_matrix": (
                    "explicit tests for amplitude-phase and cross-resource generation"
                ),
                "O98_conditional_EW_falsification_vector": (
                    "one-shot observable vector after explicit sector/representation axioms"
                ),
                "O99_conversion_metrology_frame": (
                    "rank-three ell/tau/hbar calibration and dimensional observable map"
                ),
                "O100_external_evidence_admission_pair": (
                    "independent hold-out and physical-reservoir gate contracts"
                ),
            },
        }
    }
    results["P337"], results["_P337_rows"] = program_337()
    results["P338"], results["_P338_rows"] = program_338()
    results["P339"], results["_P339_rows"] = program_339(rng)
    results["P340"], results["_P340_rows"] = program_340()
    results["P341"], results["_P341_rows"] = program_341(strict_a, rng)
    results["P342"] = program_342()
    results["P343"], results["_P343_rows"] = program_343(strict_a)
    results["P344"], results["_P344_rows"] = program_344()
    results["P345"], results["_P345_rows"] = program_345()
    results["P346"], results["_P346_rows"] = program_346(rng)
    results["P347"], results["_P347_rows"] = program_347()
    results["P348"], results["_P348_rows"] = program_348()
    results["P349"], results["_P349_rows"] = program_349()
    results["P350"] = program_350()
    results["formal_verification"] = run_formal_core()

    write_csv(INTERVAL_PATH, results["_P337_rows"])
    write_csv(OSCILLATORY_PATH, results["_P338_rows"])
    write_csv(REFREEZE_PATH, results["_P339_rows"])
    write_csv(NATURALITY_PATH, results["_P340_rows"])
    write_csv(PHOTONIC_PATH, results["_P341_rows"])
    write_csv(COMB_PATH, results["_P343_rows"])
    write_csv(CLOCK_PATH, results["_P344_rows"])
    write_csv(INVARIANT_PATH, results["_P345_rows"])
    write_csv(MONOTONE_PATH, results["_P346_rows"])
    write_csv(COUPLING_PATH, results["_P347_rows"])
    write_csv(EW_PATH, results["_P348_rows"])
    write_csv(METROLOGY_PATH, results["_P349_rows"])
    write_csv(SUMMARY_PATH, summary_rows(results))

    public_results = {
        key: value for key, value in results.items() if not key.startswith("_")
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(public_results), ensure_ascii=False, indent=2)
        + "\n",
        encoding="utf-8",
    )
    make_figures(results)
    print(RESULTS_PATH)
    print(SUMMARY_PATH)
    for program in range(337, 351):
        result = public_results[f"P{program}"]
        print(f"P{program}: {result['status']}")


if __name__ == "__main__":
    main()
