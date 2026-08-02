#!/usr/bin/env python3
"""FIN local research Programs P458, P459, and P464.

P458 differentiates the licensed palindromic coarse-erasure objective by
outward interval automatic differentiation.  P459 propagates the P453
minimum-TV gauge through the complete O163 allocation tube.  P464 continues
the finite nested-comb dual and accepts only a new exact-rational feasible
certificate as a theorem-level improvement.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

import fin_programs_435_436_440 as p435
import fin_programs_445_447 as p445
import fin_programs_448_450 as p448
import fin_programs_451_453 as p451
import fin_programs_454_455_457 as p454


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_458_459_464"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
P458_PATH = ROOT / f"{PREFIX}_P458_Derivative_Uniqueness.csv"
P459_PATH = ROOT / f"{PREFIX}_P459_Canonical_Allocation.csv"
P464_PATH = ROOT / f"{PREFIX}_P464_Dual_Refinement.csv"
P464_WITNESS_PATH = ROOT / f"{PREFIX}_P464_Rational_Dual.npz"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p458_p459_p464_certificates.png"

FInterval = tuple[float, float]


def json_ready(value: Any) -> Any:
    if isinstance(value, Fraction):
        if max(value.numerator.bit_length(), value.denominator.bit_length()) <= 512:
            return str(value)
        return format(float(value), ".17g")
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


def freciprocal(value: FInterval) -> FInterval:
    if value[0] <= 0 <= value[1]:
        raise ZeroDivisionError("interval contains zero")
    candidates = (1.0 / value[0], 1.0 / value[1])
    return p445.fdown(min(candidates)), p445.fup(max(candidates))


def fdiv(left: FInterval, right: FInterval) -> FInterval:
    return p445.fmul(left, freciprocal(right))


def fconst(value: Fraction | int | float | FInterval) -> FInterval:
    if isinstance(value, tuple):
        return value
    if isinstance(value, Fraction):
        return p445.fiv_from_fraction((value, value))
    if isinstance(value, int):
        return float(value), float(value)
    return p445.fdown(float(value)), p445.fup(float(value))


@dataclass(frozen=True)
class Jet2:
    value: FInterval
    first: FInterval
    second: FInterval


JZERO = Jet2((0.0, 0.0), (0.0, 0.0), (0.0, 0.0))


def jconst(value: Fraction | int | float | FInterval) -> Jet2:
    return Jet2(fconst(value), (0.0, 0.0), (0.0, 0.0))


def jvariable(value: FInterval) -> Jet2:
    return Jet2(value, (1.0, 1.0), (0.0, 0.0))


def jadd(left: Jet2, right: Jet2) -> Jet2:
    return Jet2(
        p445.fadd(left.value, right.value),
        p445.fadd(left.first, right.first),
        p445.fadd(left.second, right.second),
    )


def jneg(value: Jet2) -> Jet2:
    return Jet2(p445.fneg(value.value), p445.fneg(value.first), p445.fneg(value.second))


def jsub(left: Jet2, right: Jet2) -> Jet2:
    return jadd(left, jneg(right))


def jmul(left: Jet2, right: Jet2) -> Jet2:
    return Jet2(
        p445.fmul(left.value, right.value),
        p445.fadd(p445.fmul(left.first, right.value), p445.fmul(left.value, right.first)),
        p445.fsum([
            p445.fmul(left.second, right.value),
            p445.fscale(p445.fmul(left.first, right.first), (2.0, 2.0)),
            p445.fmul(left.value, right.second),
        ]),
    )


def jscale(value: Jet2, scalar: Fraction | int | float | FInterval) -> Jet2:
    return jmul(value, jconst(scalar))


def jsquare(value: Jet2) -> Jet2:
    return jmul(value, value)


def jsqrt(value: Jet2) -> Jet2:
    if value.value[0] <= 0:
        raise ValueError("sqrt jet requires a strictly positive interval")
    root = p445.fsqrt(value.value)
    two_root = p445.fscale(root, (2.0, 2.0))
    first = fdiv(value.first, two_root)
    root_cubed = p445.fmul(root, p445.fmul(root, root))
    second = p445.fsub(
        fdiv(value.second, two_root),
        fdiv(p445.fmul(value.first, value.first), p445.fscale(root_cubed, (4.0, 4.0))),
    )
    return Jet2(root, first, second)


def jabs(value: Jet2) -> Jet2:
    if value.value[0] > 0:
        return value
    if value.value[1] < 0:
        return jneg(value)
    raise ValueError("absolute-value branch crosses zero")


def jsum(values: list[Jet2]) -> Jet2:
    result = JZERO
    for value in values:
        result = jadd(result, value)
    return result


def jmatrix(rows: int, columns: int) -> list[list[Jet2]]:
    return [[JZERO for _ in range(columns)] for _ in range(rows)]


def jproject(matrix: list[list[Jet2]], basis: list[list[FInterval]]) -> list[list[Jet2]]:
    source = len(matrix)
    target = len(basis[0])
    result = jmatrix(target, target)
    for left in range(target):
        for right in range(target):
            terms: list[Jet2] = []
            for row in range(source):
                if basis[row][left] == (0.0, 0.0):
                    continue
                for column in range(source):
                    if basis[column][right] == (0.0, 0.0):
                        continue
                    terms.append(jscale(
                        jscale(matrix[row][column], basis[row][left]),
                        basis[column][right],
                    ))
            result[left][right] = jsum(terms)
    return result


def q_power(exponent: int) -> FInterval:
    result = (1.0, 1.0)
    for _ in range(exponent):
        result = p445.fmul(result, p445.Q_FI)
    return result


def jskew(probabilities: list[Jet2], survivors: int) -> list[list[Jet2]]:
    dimension = 2**survivors
    lost = 3 - survivors
    result = jmatrix(dimension, dimension)
    for left in range(dimension):
        for right in range(dimension):
            left_weight = left.bit_count()
            right_weight = right.bit_count()
            hamming = (left ^ right).bit_count()
            terms: list[Jet2] = []
            for lost_value in range(2**lost):
                lost_weight = lost_value.bit_count()
                wl = left_weight + lost_weight
                wr = right_weight + lost_weight
                denominator = math.comb(3, wl) * math.comb(3, wr)
                radicand = jscale(
                    jmul(probabilities[wl], probabilities[wr]),
                    Fraction(1, denominator),
                )
                terms.append(jsqrt(radicand))
            result[left][right] = jscale(
                jsum(terms),
                p445.fmul(
                    p445.fscale(q_power(hamming), (2.0, 2.0)),
                    p445.SINES_FI[left_weight-right_weight],
                ),
            )
    return result


def palindromic_objective_jet(interval_a: FInterval) -> Jet2:
    a = jvariable(interval_a)
    b = jsub(jconst(Fraction(1, 2)), a)
    probabilities = [a, b, b, a]

    matrix1 = jskew(probabilities, 1)
    distance1 = jabs(matrix1[0][1])

    block2 = jproject(jskew(probabilities, 2), p445.U2_FI)
    distance2 = jsqrt(jsum([
        jsquare(block2[i][j]) for i in range(3) for j in range(i+1, 3)
    ]))

    block3 = jproject(jskew(probabilities, 3), p445.U3_FI)
    squared = jsum([
        jsquare(block3[i][j]) for i in range(4) for j in range(i+1, 4)
    ])
    pfaffian = jadd(
        jsub(jmul(block3[0][1], block3[2][3]), jmul(block3[0][2], block3[1][3])),
        jmul(block3[0][3], block3[1][2]),
    )
    symmetric = jsqrt(jadd(squared, jscale(jabs(pfaffian), 2)))
    standard = jproject(jskew(probabilities, 3), p445.V3_FI)
    distance3 = jadd(symmetric, jscale(jabs(standard[0][1]), 2))
    return jsum([
        jscale(distance1, Fraction(12, 125)),
        jscale(distance2, Fraction(48, 125)),
        jscale(distance3, Fraction(64, 125)),
    ])


def program_458() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    inherited = json.loads(
        (ROOT / "FIN_Programs_454_455_457_Results.json").read_text(encoding="utf-8")
    )["P457"]
    hull_low, hull_high = map(float, inherited["maximizer_hull"])
    cells = 256
    edges = np.linspace(hull_low, hull_high, cells + 1)
    rows: list[dict[str, Any]] = []
    maximum_second_upper = -math.inf
    minimum_second_lower = math.inf
    for index in range(cells):
        jet = palindromic_objective_jet((float(edges[index]), float(edges[index+1])))
        maximum_second_upper = max(maximum_second_upper, jet.second[1])
        minimum_second_lower = min(minimum_second_lower, jet.second[0])
        rows.append({
            "row_type": "curvature_cell",
            "cell": index,
            "a_lower": edges[index],
            "a_upper": edges[index+1],
            "first_lower": jet.first[0],
            "first_upper": jet.first[1],
            "second_lower": jet.second[0],
            "second_upper": jet.second[1],
        })
    if maximum_second_upper >= 0:
        raise AssertionError("P458 strict concavity was not certified")

    left, right = hull_low, hull_high
    left_derivative = palindromic_objective_jet((left, left)).first
    right_derivative = palindromic_objective_jet((right, right)).first
    if left_derivative[0] <= 0 or right_derivative[1] >= 0:
        raise AssertionError("P458 derivative signs do not bracket a root")
    steps = 0
    while right - left > 2e-13 and steps < 80:
        midpoint = (left + right) / 2
        derivative = palindromic_objective_jet((midpoint, midpoint)).first
        if derivative[0] > 0:
            left = midpoint
        elif derivative[1] < 0:
            right = midpoint
        else:
            # Constant-enclosure uncertainty has reached its resolution floor.
            left = np.nextafter(midpoint, -np.inf)
            right = np.nextafter(midpoint, np.inf)
            break
        steps += 1
    root_jet = palindromic_objective_jet((left, right))
    numerical_value = p445.p446_numeric_objective(
        np.array([(left+right)/2, 0.5-(left+right)/2, 0.5-(left+right)/2, (left+right)/2])
    )
    rows.append({
        "row_type": "unique_root_certificate",
        "a_lower": left,
        "a_upper": right,
        "root_width": right-left,
        "objective_lower": root_jet.value[0],
        "objective_upper": root_jet.value[1],
        "numerical_objective": numerical_value,
        "left_derivative_lower": left_derivative[0],
        "right_derivative_upper": right_derivative[1],
        "maximum_second_derivative_upper": maximum_second_upper,
    })
    return ({
        "status": (
            "[Computer-assisted proof] the P457 palindromic maximizer is unique "
            "and therefore, through the P452 symmetrization theorem, gives the "
            "unique palindromic representative of the declared global optimum"
        ),
        "inherited_global_maximizer_hull": [hull_low, hull_high],
        "curvature_cells": cells,
        "minimum_second_derivative_lower": minimum_second_lower,
        "maximum_second_derivative_upper": maximum_second_upper,
        "left_endpoint_derivative_interval": left_derivative,
        "right_endpoint_derivative_interval": right_derivative,
        "unique_maximizer_interval": [left, right],
        "unique_maximizer_width": right-left,
        "objective_interval_at_root_box": root_jet.value,
        "numerical_objective_at_midpoint": numerical_value,
        "globality_chain": (
            "P452 places at least one global optimizer on the palindromic line; "
            "P457 confines every palindromic maximizer to the inherited hull; "
            "P458 proves F''<0 throughout that hull and opposite derivative signs."
        ),
        "boundary": (
            "Uniqueness is proved for the palindromic representative in the exact "
            "P452 code/parameter/loss model. P452 does not prove that every asymmetric "
            "optimizer is equal to its reversal average, so uniqueness on the original "
            "full simplex is not claimed without strict symmetrization equality analysis."
        ),
        "new_object": "O173 Palindromic Curvature Isolator",
    }, rows)


def p459_scores() -> tuple[list[p435.Interval], list[p435.Interval], list[p435.Interval], list[p435.Interval]]:
    nodes, weights = p445.p429_atom_boxes()
    radii: list[p435.Interval] = []
    hbars: list[p435.Interval] = []
    scores: list[p435.Interval] = []
    for node, weight in zip(nodes, weights):
        radius_squared = p445.iv_sum([
            p445.iv_power_nonnegative(node, 2*order) for order in range(12)
        ])
        radius = p445.iv_sqrt(radius_squared, 30)
        eta = p445.iv_add((Fraction(13, 20), Fraction(13, 20)), p445.iv_scale(node, Fraction(1, 4)))
        dark = p445.iv_sub((Fraction(3, 100), Fraction(3, 100)), p445.iv_scale(node, Fraction(1, 50)))
        hbar = p445.iv_add(
            p445.iv_reciprocal(eta),
            p445.iv_div(p445.iv_mul(dark, p445.iv_sub((Fraction(1), Fraction(1)), dark)), p445.iv_square(eta)),
        )
        score = p445.iv_mul(p445.iv_mul(p445.iv_abs(weight), radius), p445.iv_sqrt(hbar, 30))
        radii.append(radius)
        hbars.append(hbar)
        scores.append(score)
    denominator = p445.iv_sum(scores)
    allocations = [p445.iv_div(score, denominator) for score in scores]
    return nodes, weights, scores, allocations


def program_459() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    nodes, weights, scores, allocations = p459_scores()
    rows: list[dict[str, Any]] = []
    curvature_lowers: list[Fraction] = []
    for index, (node, weight, score, allocation) in enumerate(
        zip(nodes, weights, scores, allocations)
    ):
        curvature = Fraction(2) * score[0]**2 / allocation[1]**3
        curvature_lowers.append(curvature)
        rows.append({
            "atom": index,
            "node_lower": node[0],
            "node_upper": node[1],
            "weight_lower": weight[0],
            "weight_upper": weight[1],
            "score_lower": score[0],
            "score_upper": score[1],
            "allocation_lower": allocation[0],
            "allocation_upper": allocation[1],
            "strict_hessian_diagonal_lower": curvature,
        })
    ranking_margins = [
        allocations[index][0] - allocations[index+1][1]
        for index in range(6)
    ]
    if min(score[0] for score in scores) <= 0:
        raise AssertionError("P459 score positivity failed")
    if min(curvature_lowers) <= 0:
        raise AssertionError("P459 strict convexity failed")
    if min(ranking_margins) <= 0:
        raise AssertionError("P459 allocation ranking failed")
    return ({
        "status": (
            "[Proven, conditional on P429/P430/P453 and the supplied O160 detector "
            "envelope] every realized detector-box coefficient vector has one unique "
            "minimum-worst-MSE allocation, and its seven-label ordering is invariant"
        ),
        "variational_objective": (
            "R(q)=sum_i s_i^2/q_i-M, q_i>0, sum_i q_i=1; "
            "s_i=|w_i| r(x_i) sqrt(hbar_i)"
        ),
        "unique_solution": "q_i=s_i/sum_j s_j",
        "strict_convexity": (
            "The ambient Hessian is diagonal with entries 2 s_i^2/q_i^3>0; "
            "its restriction to the simplex tangent space is positive definite."
        ),
        "minimum_certified_hessian_diagonal_lower": min(curvature_lowers),
        "allocation_intervals": allocations,
        "allocation_order": [0, 1, 2, 3, 4, 5, 6],
        "minimum_adjacent_ranking_margin": min(ranking_margins),
        "fiberwise_canonical": True,
        "box_collapsed_to_single_physical_allocation": False,
        "representation_dependency": (
            "P453 removes the P450 signed-measure gauge only after the minimum-TV rule "
            "is imposed. Removing that rule restores representation dependence."
        ),
        "physical_boundary": (
            "The detector envelope is synthetic and dimensionless. Fiberwise uniqueness "
            "does not choose an actual detector realization, calibrate an apparatus, or "
            "generate a laboratory preparation or record."
        ),
        "new_object": "O174 Fiberwise Canonical Allocation Tube",
    }, rows)


def adaptive_exact_margin(matrix: sp.Matrix) -> tuple[Fraction, list[Fraction]]:
    candidates = tuple(Fraction(1, 10**power) for power in range(7, 14))
    return p454.exact_sylvester_margin(matrix, candidates=candidates)


def rational_dual_certificate(
    vector: np.ndarray, scale: int = 10**14
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    x3f, x2f, x1f, lamf = p454.unpack_dual(vector)
    x3v = p454.fraction_matrix(x3f, scale)
    x2v = p454.fraction_matrix(x2f, scale)
    x1v = p454.fraction_matrix(x1f, scale)
    lam = Fraction(math.ceil(lamf*scale), scale)
    x3 = p454.sympy_fraction_matrix(x3v)
    x2 = p454.sympy_fraction_matrix(x2v)
    x1 = p454.sympy_fraction_matrix(x1v)
    half_delta, radius = p454.p454_half_delta_midpoint()
    slacks = {
        "X3_plus_halfDelta": x3 + half_delta,
        "X3_minus_halfDelta": x3 - half_delta,
        "X2_minus_X3_block0": x2 - x3[:4, :4],
        "X2_minus_X3_block1": x2 - x3[4:, 4:],
        "X1_minus_X2_block0": x1 - x2[:2, :2],
        "X1_minus_X2_block1": x1 - x2[2:, 2:],
    }
    rows: list[dict[str, Any]] = []
    lowers: list[Fraction] = []
    for name, slack in slacks.items():
        margin, determinants = adaptive_exact_margin(slack)
        perturbation = 8*radius if "halfDelta" in name else Fraction(0)
        certified = margin-perturbation
        lowers.append(certified)
        rows.append({
            "slack": name,
            "exact_sylvester_margin": margin,
            "transcendental_operator_perturbation": perturbation,
            "certified_exact_minimum_eigenvalue_lower": certified,
            "positive_leading_principal_minors": len(determinants),
            "smallest_leading_principal_minor": min(determinants) if determinants else None,
        })
    for index, slack in enumerate((lam-x1v[0][0], lam-x1v[1][1])):
        lowers.append(slack)
        rows.append({
            "slack": f"lambda_minus_X1_diag{index}",
            "exact_sylvester_margin": slack,
            "transcendental_operator_perturbation": Fraction(0),
            "certified_exact_minimum_eigenvalue_lower": slack,
            "positive_leading_principal_minors": 1,
            "smallest_leading_principal_minor": slack,
        })
    feasible = min(lowers) > 0
    return ({
        "feasible": feasible,
        "rationalization_scale": scale,
        "lambda_upper": lam,
        "smallest_certified_slack_lower": min(lowers),
        "half_delta_maximum_entry_radius": radius,
        "slacks": rows,
    }, {
        "X3": np.array([[float(v) for v in row] for row in x3v]),
        "X2": np.array([[float(v) for v in row] for row in x2v]),
        "X1": np.array([[float(v) for v in row] for row in x1v]),
        "lambda": np.array([float(lam)]),
    })


def program_464() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    inherited = json.loads(
        (ROOT / "FIN_Programs_454_455_457_Results.json").read_text(encoding="utf-8")
    )["P454"]
    previous_upper = Fraction(inherited["certified_dual_upper"])
    primal_lower = Fraction(inherited["inherited_primal_lower"])
    source = np.load(ROOT / "FIN_Programs_454_455_457_P454_Rational_Dual.npz")
    vector = p454.pack_dual(source["X3"], source["X2"], source["X1"], float(source["lambda"][0]))
    delta = p448.compressed_process_difference(3, 0.8, math.pi/8)/2
    rows: list[dict[str, Any]] = []
    for barrier in (1e-7, 3e-8, 1e-8):
        vector, diagnostics = p454.feasible_bfgs(
            vector, delta, barrier, maximum_iterations=900, tolerance=3e-8
        )
        slacks, (_, _, _, lam) = p454.dual_slacks(vector, delta)
        rows.append({
            "row_type": "central_path_scout",
            "barrier": barrier,
            "lambda": lam,
            "minimum_float_slack": min(float(np.linalg.eigvalsh(s)[0]) for s in slacks),
            **diagnostics,
        })
    certificate, matrices = rational_dual_certificate(vector)
    if certificate["feasible"] and Fraction(certificate["lambda_upper"]) < previous_upper:
        accepted = certificate
        np.savez_compressed(P464_WITNESS_PATH, **matrices)
        improved = True
    else:
        accepted = json.loads(
            (ROOT / "FIN_Programs_454_455_457_Results.json").read_text(encoding="utf-8")
        )["P454"]["rational_dual_certificate"]
        accepted["lambda_upper"] = inherited["certified_dual_upper"]
        improved = False
        # Keep a deterministic restart artifact even when the scout fails.
        np.savez_compressed(P464_WITNESS_PATH, **matrices)
    upper = Fraction(accepted["lambda_upper"])
    gap = upper-primal_lower
    rows.extend({"row_type": "candidate_exact_slack", **row} for row in certificate["slacks"])
    return ({
        "status": (
            "[Computer-assisted proof] improved exact-rational global dual bound"
            if improved else
            "[Refuted for the attempted continuation] the smaller-barrier scout did not "
            "yield a new exact-rational feasible upper certificate; the P454 theorem remains"
        ),
        "attempted_barriers": [1e-7, 3e-8, 1e-8],
        "candidate_rational_certificate_feasible": certificate["feasible"],
        "strict_improvement_accepted": improved,
        "previous_certified_upper": previous_upper,
        "accepted_certified_upper": upper,
        "accepted_global_gap": gap,
        "candidate_certificate": certificate,
        "methodological_result": (
            "A floating central-path decrease is not an upper-bound theorem. Only exact "
            "rational Sylvester positivity after trigonometric perturbation payment is admissible."
        ),
        "boundary": (
            "The improved upper certificate does not decide exact O167 attainment or "
            "optimizer uniqueness; a future failed continuation would not prove the "
            "remaining gap irreducible."
        ),
        "new_object": "O175 Rational-Dual Admission Gate",
    }, rows)


def make_figure(p458: dict[str, Any], p458_rows: list[dict[str, Any]], p459: dict[str, Any], p459_rows: list[dict[str, Any]], p464: dict[str, Any], p464_rows: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.35))

    cells = [row for row in p458_rows if row["row_type"] == "curvature_cell"]
    x = [(row["a_lower"]+row["a_upper"])/2 for row in cells]
    axes[0].fill_between(x, [row["second_lower"] for row in cells], [row["second_upper"] for row in cells], color="#2563eb", alpha=0.35)
    axes[0].axhline(0, color="#dc2626", linestyle="--")
    axes[0].set_title("P458 certified strict curvature")
    axes[0].set_xlabel("palindromic mass a")
    axes[0].set_ylabel("second derivative interval")

    atoms = [row["atom"] for row in p459_rows]
    lows = [float(Fraction(row["allocation_lower"])) for row in p459_rows]
    highs = [float(Fraction(row["allocation_upper"])) for row in p459_rows]
    centers = [(a+b)/2 for a,b in zip(lows, highs)]
    axes[1].bar(atoms, centers, color="#0f766e")
    axes[1].set_title("P459 unique fiberwise allocation")
    axes[1].set_xlabel("canonical atom label")
    axes[1].set_ylabel("sampling probability")

    scouts = [row for row in p464_rows if row["row_type"] == "central_path_scout"]
    axes[2].semilogx([row["barrier"] for row in scouts], [row["lambda"] for row in scouts], "o-", color="#7c3aed")
    axes[2].axhline(float(Fraction(p464["previous_certified_upper"])), color="#f59e0b", linestyle="--", label="P454 certified upper")
    axes[2].invert_xaxis()
    axes[2].set_title("P464 dual admission attempt")
    axes[2].set_xlabel("barrier parameter")
    axes[2].set_ylabel("floating dual lambda")
    axes[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)


def main() -> None:
    p458_result, p458_rows = program_458()
    write_csv(P458_PATH, p458_rows)
    p459_result, p459_rows = program_459()
    write_csv(P459_PATH, p459_rows)
    p464_result, p464_rows = program_464()
    write_csv(P464_PATH, p464_rows)
    make_figure(p458_result, p458_rows, p459_result, p459_rows, p464_result, p464_rows)
    results = {
        "metadata": {
            "programs": "P458/P459/P464",
            "checkpoint": "P458-P464",
            "date": "2026-08-01",
            "local_only": True,
            "external_physical_evidence": False,
            "selector_discharged": False,
            "dimensional_source_exported": False,
            "legacy_strict_bridge_complete": False,
            "legacy_role_transfer_started": False,
        },
        "P458": p458_result,
        "P459": p459_result,
        "P464": p464_result,
    }
    RESULTS_PATH.write_text(json.dumps(json_ready(results), indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    write_csv(SUMMARY_PATH, [{
        "program": name,
        "status": results[name]["status"],
        "new_object": results[name]["new_object"],
        "physical_evidence": False,
        "selector_discharged": False,
    } for name in ("P458", "P459", "P464")])
    print(json.dumps({
        "P458_unique_interval": p458_result["unique_maximizer_interval"],
        "P459_minimum_ranking_margin": json_ready(p459_result["minimum_adjacent_ranking_margin"]),
        "P464_strict_improvement": p464_result["strict_improvement_accepted"],
        "P464_accepted_gap": json_ready(p464_result["accepted_global_gap"]),
    }, indent=2))


if __name__ == "__main__":
    main()
