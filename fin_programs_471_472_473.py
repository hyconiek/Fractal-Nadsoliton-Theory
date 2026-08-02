#!/usr/bin/env python3
"""FIN local research Programs P471--P473.

The batch replaces the floating matrix-functional KKT locator by a finite
polynomial Riccati system, constructs an exact rational interval evaluator,
and attempts a Krawczyk proof of exact O167 primal--dual attainment.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
import sympy as sp

import fin_programs_448_450 as p448
import fin_programs_454_455_457 as p454
import fin_programs_465_468_469 as p465


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_471_472_473"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
P471_PATH = ROOT / f"{PREFIX}_P471_Polynomial_Core.csv"
P472_PATH = ROOT / f"{PREFIX}_P472_Interval_Provider.csv"
P473_PATH = ROOT / f"{PREFIX}_P473_Krawczyk.csv"
P473_WITNESS_PATH = ROOT / f"{PREFIX}_P473_Root_Box.npz"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p471_p473_polynomial_krawczyk.png"

Interval = tuple[Fraction, Fraction]


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


def polynomial_system() -> dict[str, Any]:
    """Construct the exact 13-variable Riccati/KKT polynomial system."""

    A, B, u = sp.symbols("A B u", real=True)
    L, a, b, c, d, e, f, g, h, i = sp.symbols(
        "L a b c d e f g h i", real=True
    )
    s1, s2, s3 = sp.symbols("s1 s2 s3", real=True)
    variables = (A, B, u, L, a, b, c, d, e, f, g, h, i)
    sine_symbols = (s1, s2, s3)

    C = sp.Rational(1, 2) - A - 2 * B
    left = sp.Matrix([
        [A, u, u, 0],
        [u, B, 0, -u],
        [u, 0, B, -u],
        [0, -u, -u, C],
    ])
    right = sp.Matrix([
        [C, -u, -u, 0],
        [-u, B, 0, u],
        [-u, 0, B, u],
        [0, u, u, A],
    ])
    normalizer = sp.diag(left, right)

    x2 = sp.Matrix([
        [L, a, a, c],
        [a, L, b, a],
        [a, b, L, a],
        [c, a, a, L],
    ])
    cross = sp.Matrix([
        [d, e, e, f],
        [g, h, h, e],
        [g, h, h, e],
        [i, g, g, d],
    ])
    x3 = x2.row_join(cross).col_join(cross.T.row_join(x2))

    q = sp.Rational(4, 5)
    sine = {0: 0, 1: s1, 2: s2, 3: s3, -1: -s1, -2: -s2, -3: -s3}
    delta = sp.zeros(8)
    for row_index in range(8):
        for column_index in range(8):
            exponent = (row_index ^ column_index).bit_count()
            difference = row_index.bit_count() - column_index.bit_count()
            delta[row_index, column_index] = (
                2 * sp.I * q**exponent * sine[difference]
            )

    residual = sp.expand(x3 * normalizer * x3 - delta * normalizer * delta / 4)
    groups: list[tuple[sp.Expr, list[tuple[int, int]]]] = []
    for row_index in range(8):
        for column_index in range(row_index, 8):
            expression = sp.expand(residual[row_index, column_index])
            for group_index, (representative, positions) in enumerate(groups):
                if sp.expand(expression - representative) == 0:
                    positions.append((row_index, column_index))
                    groups[group_index] = (representative, positions)
                    break
            else:
                groups.append((expression, [(row_index, column_index)]))
    if len(groups) != 13:
        raise AssertionError(f"expected 13 residual orbits, found {len(groups)}")
    equations = tuple(group[0] for group in groups)
    jacobian = sp.Matrix(equations).jacobian(variables)
    all_symbols = variables + sine_symbols
    equation_polys = tuple(sp.Poly(value, *all_symbols, domain=sp.QQ) for value in equations)
    jacobian_polys = tuple(
        tuple(sp.Poly(jacobian[row_index, column_index], *all_symbols, domain=sp.QQ)
              for column_index in range(13))
        for row_index in range(13)
    )
    return {
        "variables": variables,
        "sines": sine_symbols,
        "normalizer": normalizer,
        "x3": x3,
        "delta": delta,
        "equations": equations,
        "jacobian": jacobian,
        "equation_polys": equation_polys,
        "jacobian_polys": jacobian_polys,
        "groups": [positions for _, positions in groups],
    }


def sine_midpoints() -> tuple[float, float, float]:
    return tuple(
        math.sin(math.pi * index / 8) for index in (1, 2, 3)
    )


def initial_polynomial_vector() -> np.ndarray:
    p469 = json.loads(
        (ROOT / "FIN_Programs_465_468_469_Results.json").read_text(encoding="utf-8")
    )["P469"]
    parameters = np.array(p469["stationary_parameters"], dtype=float)
    delta = p448.compressed_process_difference(3, 0.8, math.pi / 8)
    _, x3, _, _, _, _, _ = p465.o167_support_ladder(parameters, delta)
    x3 = np.real(x3)
    support_parameters = [
        x3[0, 0], x3[0, 1], x3[1, 2], x3[0, 3],
        x3[0, 4], x3[0, 5], x3[0, 7],
        x3[1, 4], x3[1, 5], x3[3, 4],
    ]
    return np.concatenate((parameters, support_parameters))


def numerical_root(system: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arguments = system["variables"] + system["sines"]
    function = sp.lambdify(arguments, system["equations"], modules="numpy")
    jacobian_function = sp.lambdify(arguments, system["jacobian"], modules="numpy")
    sines = sine_midpoints()
    evaluate = lambda value: np.asarray(function(*value, *sines), dtype=float).reshape(13)
    evaluate_jacobian = lambda value: np.asarray(
        jacobian_function(*value, *sines), dtype=float
    ).reshape(13, 13)
    solution = root(
        evaluate,
        initial_polynomial_vector(),
        jac=evaluate_jacobian,
        method="lm",
        options={"ftol": 1e-14, "xtol": 1e-14, "gtol": 1e-14, "maxiter": 3000},
    )
    point = np.asarray(solution.x, dtype=float)
    return point, evaluate(point), evaluate_jacobian(point)


def iv_add(left: Interval, right: Interval) -> Interval:
    return left[0] + right[0], left[1] + right[1]


def iv_neg(value: Interval) -> Interval:
    return -value[1], -value[0]


def iv_sub(left: Interval, right: Interval) -> Interval:
    return iv_add(left, iv_neg(right))


def iv_mul(left: Interval, right: Interval) -> Interval:
    products = (
        left[0] * right[0], left[0] * right[1],
        left[1] * right[0], left[1] * right[1],
    )
    return min(products), max(products)


def iv_scale(value: Interval, scalar: Fraction) -> Interval:
    return iv_mul(value, (scalar, scalar))


def iv_pow(value: Interval, exponent: int) -> Interval:
    if exponent == 0:
        return Fraction(1), Fraction(1)
    result = value
    for _ in range(exponent - 1):
        result = iv_mul(result, value)
    return result


def rational_enclosure(value: Interval, scale: int) -> Interval:
    lower_scaled = value[0] * scale
    upper_scaled = value[1] * scale
    lower = lower_scaled.numerator // lower_scaled.denominator
    upper = -((-upper_scaled.numerator) // upper_scaled.denominator)
    return Fraction(lower, scale), Fraction(upper, scale)


def polynomial_interval(poly: sp.Poly, values: tuple[Interval, ...]) -> Interval:
    result = (Fraction(0), Fraction(0))
    for monomial, coefficient in poly.terms():
        rational = Fraction(int(coefficient.p), int(coefficient.q))
        term = (rational, rational)
        for value, exponent in zip(values, monomial):
            if exponent:
                term = iv_mul(term, iv_pow(value, exponent))
        result = iv_add(result, term)
    return result


def exact_sine_boxes(scale: int = 10**30) -> tuple[Interval, Interval, Interval]:
    return tuple(
        rational_enclosure(p448.p449_sine_interval(index), scale)
        for index in (1, 2, 3)
    )


def rational_center(point: np.ndarray, scale: int = 10**16) -> tuple[Fraction, ...]:
    return tuple(Fraction(round(float(value) * scale), scale) for value in point)


def rational_matrix(values: np.ndarray, scale: int = 10**14) -> list[list[Fraction]]:
    return [
        [Fraction(round(float(value) * scale), scale) for value in row]
        for row in values
    ]


def interval_matrix_product(
    left: list[list[Fraction]], right: list[list[Interval]]
) -> list[list[Interval]]:
    rows = len(left)
    shared = len(left[0])
    columns = len(right[0])
    result = [[(Fraction(0), Fraction(0)) for _ in range(columns)] for _ in range(rows)]
    for row_index in range(rows):
        for column_index in range(columns):
            cell = (Fraction(0), Fraction(0))
            for index in range(shared):
                cell = iv_add(cell, iv_scale(right[index][column_index], left[row_index][index]))
            result[row_index][column_index] = cell
    return result


def krawczyk_certificate(
    system: dict[str, Any],
    point: np.ndarray,
    jacobian_numeric: np.ndarray,
    radius: Fraction,
    center_scale: int = 10**16,
    inverse_scale: int = 10**14,
) -> dict[str, Any]:
    center = rational_center(point, center_scale)
    box = tuple((value - radius, value + radius) for value in center)
    sines = exact_sine_boxes()
    point_values = tuple((value, value) for value in center) + sines
    box_values = box + sines
    f_point = [
        polynomial_interval(poly, point_values) for poly in system["equation_polys"]
    ]
    j_box = [
        [polynomial_interval(poly, box_values) for poly in row]
        for row in system["jacobian_polys"]
    ]
    inverse = rational_matrix(np.linalg.inv(jacobian_numeric), inverse_scale)
    c_f = []
    for row_index in range(13):
        value = (Fraction(0), Fraction(0))
        for column_index in range(13):
            value = iv_add(
                value,
                iv_scale(f_point[column_index], inverse[row_index][column_index]),
            )
        c_f.append(value)
    c_j = interval_matrix_product(inverse, j_box)
    iteration = []
    inclusion_margins = []
    for row_index in range(13):
        base = iv_sub((center[row_index], center[row_index]), c_f[row_index])
        correction = (Fraction(0), Fraction(0))
        for column_index in range(13):
            identity = Fraction(1) if row_index == column_index else Fraction(0)
            coefficient = iv_sub((identity, identity), c_j[row_index][column_index])
            correction = iv_add(correction, iv_mul(coefficient, (-radius, radius)))
        image = iv_add(base, correction)
        iteration.append(image)
        inclusion_margins.append(min(
            image[0] - box[row_index][0],
            box[row_index][1] - image[1],
        ))
    included = min(inclusion_margins) > 0
    preconditioned_radius = max(
        max(abs(float(value[0])), abs(float(value[1])))
        for row in c_j
        for value in row
    )
    return {
        "included": included,
        "center": center,
        "radius": radius,
        "box": box,
        "image": iteration,
        "minimum_inclusion_margin": min(inclusion_margins),
        "maximum_image_width": max(value[1] - value[0] for value in iteration),
        "maximum_point_residual_radius": max(
            max(abs(value[0]), abs(value[1])) for value in f_point
        ),
        "preconditioner_scale": inverse_scale,
        "center_scale": center_scale,
        "sine_boxes": sines,
        "maximum_absolute_CJ_entry": preconditioned_radius,
    }


def structured_numeric_matrices(vector: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    A, B, u, L, a, b, c, d, e, f, g, h, i = vector
    C = 0.5 - A - 2 * B
    left = np.array([
        [A, u, u, 0], [u, B, 0, -u],
        [u, 0, B, -u], [0, -u, -u, C],
    ])
    right = np.array([
        [C, -u, -u, 0], [-u, B, 0, u],
        [-u, 0, B, u], [0, u, u, A],
    ])
    normalizer = np.zeros((8, 8))
    normalizer[:4, :4] = left
    normalizer[4:, 4:] = right
    x2 = np.array([
        [L, a, a, c], [a, L, b, a],
        [a, b, L, a], [c, a, a, L],
    ])
    cross = np.array([
        [d, e, e, f], [g, h, h, e],
        [g, h, h, e], [i, g, g, d],
    ])
    x3 = np.block([[x2, cross], [cross.T, x2]])
    return normalizer, x3


def program_471(
    system: dict[str, Any], point: np.ndarray, residual: np.ndarray, jacobian: np.ndarray
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    determinant = float(np.linalg.det(jacobian))
    rows = []
    for index, positions in enumerate(system["groups"]):
        rows.append({
            "row_type": "exact_residual_orbit",
            "orbit": index,
            "representative": positions[0],
            "positions": positions,
            "multiplicity": len(positions),
        })
    rows.append({
        "row_type": "floating_root",
        "point": point,
        "residual_infinity_norm": float(np.linalg.norm(residual, np.inf)),
        "jacobian_determinant": determinant,
        "jacobian_singular_values": singular_values,
    })
    return ({
        "status": (
            "[Proven] the interior O167 KKT support condition is equivalent to a "
            "13-variable polynomial Riccati system with exactly 13 residual orbits; "
            "[Strong evidence] its located root has a well-conditioned Jacobian"
        ),
        "unknowns": [str(value) for value in system["variables"]],
        "upper_triangle_equations": 36,
        "exact_residual_orbits": len(system["groups"]),
        "residual_orbits": system["groups"],
        "polynomial_degree_maximum": max(poly.total_degree() for poly in system["equation_polys"]),
        "floating_root": point,
        "floating_residual_infinity_norm": float(np.linalg.norm(residual, np.inf)),
        "jacobian_determinant": determinant,
        "jacobian_singular_values": singular_values,
        "jacobian_condition_number": float(singular_values[0] / singular_values[-1]),
        "equivalence_theorem": (
            "For N>0 and X3>0, X3*N*X3=Delta*N*Delta/4 iff "
            "N^(1/2) X3 N^(1/2)=|N^(1/2) Delta N^(1/2)|/2. "
            "The imposed 10-parameter X3 pattern is exactly the full nested "
            "block-equality KKT ladder under the declared O167 symmetries."
        ),
        "boundary": (
            "The symbolic reduction is exact, but floating residual and Jacobian "
            "conditioning alone do not prove an exact real root."
        ),
        "new_object": "O179 Polynomialized Interior KKT Core",
    }, rows)


def program_472(
    system: dict[str, Any], point: np.ndarray, jacobian: np.ndarray
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    attempts = []
    certificates = []
    for radius in (
        Fraction(1, 10**8),
        Fraction(1, 10**9),
        Fraction(1, 10**10),
        Fraction(1, 10**11),
        Fraction(1, 10**12),
        Fraction(1, 10**13),
        Fraction(3, 10**14),
        Fraction(1, 10**14),
    ):
        certificate = krawczyk_certificate(system, point, jacobian, radius)
        certificates.append(certificate)
        attempts.append({
            "row_type": "interval_provider_attempt",
            "radius": radius,
            "included": certificate["included"],
            "minimum_inclusion_margin": certificate["minimum_inclusion_margin"],
            "maximum_image_width": certificate["maximum_image_width"],
            "maximum_point_residual_radius": certificate["maximum_point_residual_radius"],
            "maximum_absolute_CJ_entry": certificate["maximum_absolute_CJ_entry"],
        })
    admitted = [value for value in certificates if value["included"]]
    selected = min(admitted, key=lambda value: value["radius"]) if admitted else certificates[-1]
    rows = list(attempts)
    for index, (box, image) in enumerate(zip(selected["box"], selected["image"])):
        rows.append({
            "row_type": "selected_interval_coordinate",
            "index": index,
            "variable": str(system["variables"][index]),
            "box_lower": box[0],
            "box_upper": box[1],
            "krawczyk_lower": image[0],
            "krawczyk_upper": image[1],
            "strict_left_margin": image[0] - box[0],
            "strict_right_margin": box[1] - image[1],
        })
    return ({
        "status": (
            "[Computer-assisted proof] exact rational interval polynomial and "
            "Jacobian provider admits at least one strict Krawczyk box"
            if admitted else
            "[Refuted for attempted boxes] the exact interval provider did not "
            "obtain a strict Krawczyk inclusion"
        ),
        "exact_rational_interval_arithmetic": True,
        "sine_enclosure_decimal_scale": 10**30,
        "tested_radii": [value["radius"] for value in certificates],
        "admitted_radii": [value["radius"] for value in admitted],
        "selected_radius": selected["radius"],
        "strict_inclusion": selected["included"],
        "minimum_inclusion_margin": selected["minimum_inclusion_margin"],
        "maximum_krawczyk_image_width": selected["maximum_image_width"],
        "maximum_point_residual_radius": selected["maximum_point_residual_radius"],
        "provider_boundary": (
            "The provider encloses the polynomial Riccati system and all Jacobian "
            "entries with exact fractions. It does not rely on interval eigendecomposition."
        ),
        "new_object": "O180 Rational Polynomial Krawczyk Provider",
    }, rows, certificates)


def rational_matrix_from_center(
    system: dict[str, Any], center: tuple[Fraction, ...], key: str
) -> sp.Matrix:
    substitution = {
        variable: sp.Rational(value.numerator, value.denominator)
        for variable, value in zip(system["variables"], center)
    }
    return sp.Matrix(system[key].subs(substitution))


def program_473(
    system: dict[str, Any],
    point: np.ndarray,
    provider: dict[str, Any],
    certificates: list[dict[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    admitted = [value for value in certificates if value["included"]]
    if not admitted:
        return ({
            "status": "[Open] no strict Krawczyk box was admitted",
            "exact_root_proved": False,
            "exact_O167_attainment_proved": False,
            "boundary": "P472 did not supply the required exact root certificate.",
            "new_object": "O181 Exact O167 Primal-Dual Contact Certificate",
        }, [])
    selected = min(admitted, key=lambda value: value["radius"])
    center = selected["center"]
    radius = selected["radius"]
    normalizer_center = rational_matrix_from_center(system, center, "normalizer")
    x3_center = rational_matrix_from_center(system, center, "x3")
    normalizer_margin, _ = p454.exact_sylvester_margin(
        normalizer_center,
        candidates=(Fraction(1, 100), Fraction(1, 1000), Fraction(1, 10000)),
    )
    x3_margin, _ = p454.exact_sylvester_margin(
        x3_center,
        candidates=(Fraction(1, 100), Fraction(1, 1000), Fraction(1, 10000)),
    )
    # Every entry of N moves by at most 5r under A,B,u variation; every entry
    # of structured X3 moves by at most r. The n*max-entry norm bound is paid.
    normalizer_perturbation = 8 * 5 * radius
    x3_perturbation = 8 * radius
    normalizer_positive_lower = normalizer_margin - normalizer_perturbation
    x3_positive_lower = x3_margin - x3_perturbation
    exact_root = (
        selected["included"]
        and normalizer_positive_lower > 0
        and x3_positive_lower > 0
    )
    lambda_index = 3
    lambda_box = selected["box"][lambda_index]
    np.savez_compressed(
        P473_WITNESS_PATH,
        center=np.array([float(value) for value in center]),
        radius=np.array([float(radius)]),
        lower=np.array([float(value[0]) for value in selected["box"]]),
        upper=np.array([float(value[1]) for value in selected["box"]]),
    )
    rows = [
        {
            "row_type": "positivity_certificate",
            "matrix": "N",
            "center_sylvester_margin": normalizer_margin,
            "box_operator_perturbation": normalizer_perturbation,
            "certified_box_positive_lower": normalizer_positive_lower,
        },
        {
            "row_type": "positivity_certificate",
            "matrix": "X3",
            "center_sylvester_margin": x3_margin,
            "box_operator_perturbation": x3_perturbation,
            "certified_box_positive_lower": x3_positive_lower,
        },
        {
            "row_type": "exact_value_interval",
            "lambda_lower": lambda_box[0],
            "lambda_upper": lambda_box[1],
            "width": lambda_box[1] - lambda_box[0],
        },
    ]
    return ({
        "status": (
            "[Computer-assisted proof] a unique exact positive polynomial KKT "
            "root exists in the admitted box and yields exact O167 primal-dual "
            "attainment for the declared three-slot ordered comb"
            if exact_root else
            "[Open] Krawczyk inclusion was obtained but box positivity was not certified"
        ),
        "exact_root_proved": exact_root,
        "root_unique_in_box": selected["included"],
        "exact_O167_attainment_proved": exact_root,
        "full_cone_optimizer_uniqueness_proved": False,
        "root_box_radius": radius,
        "exact_global_value_interval": lambda_box,
        "exact_global_value_interval_width": lambda_box[1] - lambda_box[0],
        "normalizer_center_sylvester_margin": normalizer_margin,
        "normalizer_box_positive_lower": normalizer_positive_lower,
        "x3_center_sylvester_margin": x3_margin,
        "x3_box_positive_lower": x3_positive_lower,
        "attainment_theorem": (
            "At the exact positive root, the Riccati identity makes X3 the "
            "trace-norm support and hence X3>=+/-Delta/2. The polynomial pattern "
            "makes every lower dual slack zero. Recursive trace identities give "
            "primal value Tr(N X3)=lambda, so weak duality is attained exactly."
        ),
        "boundary": (
            "The theorem is confined to the declared reduced three-slot channel "
            "and proves one exact optimal O167 normalizer. It does not by itself "
            "exclude other full-cone optimizers outside the certified root box."
        ),
        "new_object": "O181 Exact O167 Primal-Dual Contact Certificate",
    }, rows)


def make_figure(
    p471: dict[str, Any], p472: dict[str, Any], p472_rows: list[dict[str, Any]],
    p473: dict[str, Any],
) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.4))
    singular = p471["jacobian_singular_values"]
    axes[0].bar(range(1, 14), singular, color="#2563eb")
    axes[0].set_yscale("log")
    axes[0].set_title("P471 polynomial Jacobian")
    axes[0].set_xlabel("singular-value index")
    axes[0].set_ylabel("singular value")

    attempts = [row for row in p472_rows if row["row_type"] == "interval_provider_attempt"]
    radii = [float(Fraction(row["radius"])) for row in attempts]
    margins = [float(Fraction(row["minimum_inclusion_margin"])) for row in attempts]
    colors = ["#0f766e" if row["included"] else "#dc2626" for row in attempts]
    axes[1].bar([f"{value:.0e}" for value in radii], margins, color=colors)
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_title("P472 Krawczyk inclusion margin")
    axes[1].set_xlabel("box radius")
    axes[1].set_ylabel("strict inward margin")

    if p473.get("exact_root_proved"):
        lower, upper = [float(Fraction(value)) for value in p473["exact_global_value_interval"]]
        center = (lower + upper) / 2
        axes[2].errorbar(
            [0], [center], yerr=[[(center - lower)], [(upper - center)]],
            fmt="o", capsize=5, color="#7c3aed",
        )
        axes[2].ticklabel_format(axis="y", style="plain", useOffset=False)
        axes[2].set_xticks([0], ["exact root box"])
        axes[2].set_title("P473 certified optimum value")
        axes[2].set_ylabel("global half-distance")
    else:
        axes[2].text(0.5, 0.5, "exact root open", ha="center", va="center")
        axes[2].set_axis_off()
    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)


def main() -> None:
    system = polynomial_system()
    point, residual, jacobian = numerical_root(system)
    p471_result, p471_rows = program_471(system, point, residual, jacobian)
    write_csv(P471_PATH, p471_rows)
    p472_result, p472_rows, certificates = program_472(system, point, jacobian)
    write_csv(P472_PATH, p472_rows)
    p473_result, p473_rows = program_473(
        system, point, p472_result, certificates
    )
    write_csv(P473_PATH, p473_rows)
    make_figure(p471_result, p472_result, p472_rows, p473_result)
    results = {
        "metadata": {
            "programs": "P471/P472/P473",
            "checkpoint": "P471-P473",
            "date": "2026-08-01",
            "local_only": True,
            "external_physical_evidence": False,
            "selector_discharged": False,
            "dimensional_source_exported": False,
            "legacy_strict_bridge_complete": False,
            "legacy_role_transfer_started": False,
        },
        "P471": p471_result,
        "P472": p472_result,
        "P473": p473_result,
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(results), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_csv(SUMMARY_PATH, [
        {
            "program": "P471", "status": p471_result["status"],
            "metric": "residual orbits", "value": p471_result["exact_residual_orbits"],
            "new_object": p471_result["new_object"],
        },
        {
            "program": "P472", "status": p472_result["status"],
            "metric": "minimum inclusion margin", "value": p472_result["minimum_inclusion_margin"],
            "new_object": p472_result["new_object"],
        },
        {
            "program": "P473", "status": p473_result["status"],
            "metric": "exact root proved", "value": p473_result["exact_root_proved"],
            "new_object": p473_result["new_object"],
        },
    ])
    print(json.dumps(json_ready(results), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
