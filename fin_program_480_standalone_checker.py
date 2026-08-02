#!/usr/bin/env python3
"""Standard-library-only replay checker for FIN O180/O181.

Allowed imports are deliberately restricted to fractions, json, pathlib,
sys, and typing.  The checker recomputes the rational interval Krawczyk
inclusion and exact Sylvester positivity payments from serialized monomials.
It does not trust a precomputed Boolean verdict.
"""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parent
DEFAULT_CERTIFICATE = ROOT / "FIN_Program_480_Standalone_Certificate.json"
DEFAULT_RESULT = ROOT / "FIN_Program_480_Standalone_Check_Result.json"
Interval = tuple[Fraction, Fraction]


def iv_add(left: Interval, right: Interval) -> Interval:
    return left[0]+right[0], left[1]+right[1]


def iv_neg(value: Interval) -> Interval:
    return -value[1], -value[0]


def iv_sub(left: Interval, right: Interval) -> Interval:
    return iv_add(left, iv_neg(right))


def iv_mul(left: Interval, right: Interval) -> Interval:
    products = (
        left[0]*right[0], left[0]*right[1],
        left[1]*right[0], left[1]*right[1],
    )
    return min(products), max(products)


def iv_scale(value: Interval, scalar: Fraction) -> Interval:
    return iv_mul(value, (scalar, scalar))


def iv_pow(value: Interval, exponent: int) -> Interval:
    if exponent == 0:
        return Fraction(1), Fraction(1)
    result = Fraction(1), Fraction(1)
    base = value
    power = exponent
    while power:
        if power & 1:
            result = iv_mul(result, base)
        base = iv_mul(base, base)
        power //= 2
    return result


def parse_interval(value: list[str]) -> Interval:
    return Fraction(value[0]), Fraction(value[1])


def polynomial_interval(
    terms: list[dict[str, Any]], values: tuple[Interval, ...]
) -> Interval:
    total: Interval = (Fraction(0), Fraction(0))
    for term in terms:
        monomial: Interval = (Fraction(term["coefficient"]), Fraction(term["coefficient"]))
        exponents = term["exponents"]
        if len(exponents) != len(values):
            raise ValueError("polynomial exponent length mismatch")
        for value, exponent in zip(values, exponents):
            if exponent:
                monomial = iv_mul(monomial, iv_pow(value, int(exponent)))
        total = iv_add(total, monomial)
    return total


def rational_matrix(value: list[list[str]]) -> list[list[Fraction]]:
    return [[Fraction(entry) for entry in row] for row in value]


def precondition(
    left: list[list[Fraction]], right: list[list[Interval]]
) -> list[list[Interval]]:
    rows = len(left)
    shared = len(left[0])
    columns = len(right[0])
    result = [[(Fraction(0), Fraction(0)) for _ in range(columns)] for _ in range(rows)]
    for row in range(rows):
        for column in range(columns):
            value: Interval = (Fraction(0), Fraction(0))
            for inner in range(shared):
                value = iv_add(value, iv_scale(right[inner][column], left[row][inner]))
            result[row][column] = value
    return result


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    """Exact fraction Gaussian elimination with row pivoting."""

    work = [row[:] for row in matrix]
    dimension = len(work)
    sign = 1
    product = Fraction(1)
    for column in range(dimension):
        pivot = next((row for row in range(column, dimension) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        diagonal = work[column][column]
        product *= diagonal
        for row in range(column+1, dimension):
            multiplier = work[row][column]/diagonal
            for inner in range(column+1, dimension):
                work[row][inner] -= multiplier*work[column][inner]
            work[row][column] = Fraction(0)
    return sign*product


def sylvester_after_shift(
    matrix: list[list[Fraction]], margin: Fraction
) -> tuple[bool, list[Fraction]]:
    shifted = [row[:] for row in matrix]
    for index in range(len(shifted)):
        shifted[index][index] -= margin
    minors = [
        determinant([row[:size] for row in shifted[:size]])
        for size in range(1, len(shifted)+1)
    ]
    return all(value > 0 for value in minors), minors


def replay(certificate: dict[str, Any]) -> dict[str, Any]:
    if certificate.get("format") != "FIN-P480-O180-O181-v1":
        raise ValueError("unsupported P480 certificate format")
    dimension = int(certificate["variable_count"])
    if dimension != 13 or certificate["equation_count"] != 13:
        raise ValueError("unexpected polynomial dimension")
    center = tuple(Fraction(value) for value in certificate["center"])
    radius = Fraction(certificate["radius"])
    sine_boxes = tuple(parse_interval(value) for value in certificate["sine_boxes"])
    point_values = tuple((value, value) for value in center)+sine_boxes
    variable_box = tuple((value-radius, value+radius) for value in center)
    box_values = variable_box+sine_boxes
    equations = certificate["equations"]
    jacobian = certificate["jacobian"]
    f_point = [polynomial_interval(poly, point_values) for poly in equations]
    j_box = [
        [polynomial_interval(poly, box_values) for poly in row]
        for row in jacobian
    ]
    inverse = rational_matrix(certificate["preconditioner"])
    c_f: list[Interval] = []
    for row in range(dimension):
        value: Interval = (Fraction(0), Fraction(0))
        for column in range(dimension):
            value = iv_add(value, iv_scale(f_point[column], inverse[row][column]))
        c_f.append(value)
    c_j = precondition(inverse, j_box)
    images: list[Interval] = []
    margins: list[Fraction] = []
    contraction_rows: list[Fraction] = []
    for row in range(dimension):
        base = iv_sub((center[row], center[row]), c_f[row])
        correction: Interval = (Fraction(0), Fraction(0))
        row_sum = Fraction(0)
        for column in range(dimension):
            identity = Fraction(int(row == column))
            coefficient = iv_sub((identity, identity), c_j[row][column])
            correction = iv_add(correction, iv_mul(coefficient, (-radius, radius)))
            row_sum += max(abs(coefficient[0]), abs(coefficient[1]))
        image = iv_add(base, correction)
        images.append(image)
        margins.append(min(
            image[0]-variable_box[row][0],
            variable_box[row][1]-image[1],
        ))
        contraction_rows.append(row_sum)
    krawczyk_included = min(margins) > 0

    margin = Fraction(certificate["sylvester_test_margin"])
    normalizer = rational_matrix(certificate["normalizer_center"])
    x3 = rational_matrix(certificate["X3_center"])
    normalizer_center_positive, normalizer_minors = sylvester_after_shift(normalizer, margin)
    x3_center_positive, x3_minors = sylvester_after_shift(x3, margin)
    matrix_dimension = int(certificate["matrix_dimension"])
    normalizer_perturbation = (
        matrix_dimension*int(certificate["normalizer_entry_lipschitz"])*radius
    )
    x3_perturbation = (
        matrix_dimension*int(certificate["X3_entry_lipschitz"])*radius
    )
    normalizer_lower = margin-normalizer_perturbation
    x3_lower = margin-x3_perturbation
    positivity = (
        normalizer_center_positive and x3_center_positive
        and normalizer_lower > 0 and x3_lower > 0
    )
    lambda_index = int(certificate["lambda_coordinate"])
    value_interval = variable_box[lambda_index]
    passed = krawczyk_included and positivity
    return {
        "format": certificate["format"],
        "status": (
            "[Computer-assisted proof replayed] exact O180 Krawczyk root and "
            "O181 box positivity verified with the Python standard library"
            if passed else
            "[Failed] at least one exact P480 replay condition was rejected"
        ),
        "passed": passed,
        "krawczyk_included": krawczyk_included,
        "minimum_inclusion_margin": str(min(margins)),
        "maximum_contraction_row_sum": str(max(contraction_rows)),
        "maximum_point_residual_radius": str(max(
            max(abs(value[0]), abs(value[1])) for value in f_point
        )),
        "normalizer_center_shifted_sylvester": normalizer_center_positive,
        "X3_center_shifted_sylvester": x3_center_positive,
        "normalizer_leading_minor_count": len(normalizer_minors),
        "X3_leading_minor_count": len(x3_minors),
        "normalizer_box_positive_lower": str(normalizer_lower),
        "X3_box_positive_lower": str(x3_lower),
        "exact_global_value_interval": [str(value_interval[0]), str(value_interval[1])],
        "exact_global_value_interval_width": str(value_interval[1]-value_interval[0]),
        "boundary": (
            "The checker verifies the finite exact root and positivity certificate. "
            "The Riccati-to-support and trace-telescope implications remain analytic "
            "theorems; no physical, selector, dimensional, or laboratory claim is made."
        ),
    }


def main() -> None:
    certificate_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CERTIFICATE
    result_path = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_RESULT
    certificate = json.loads(certificate_path.read_text(encoding="utf-8"))
    result = replay(certificate)
    result_path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
