#!/usr/bin/env python3
"""Independent exact-arithmetic checker for the P351/P352 certificates.

The checker deliberately imports no FIN research module.  It reimplements
the rational interval, Krawczyk, Bernstein, Taylor, fifth-root, and
Vandermonde calculations with the Python standard library, then compares the
resulting endpoint identities with the archived Release 10.31 CSV files.

This gives implementation diversity inside one language/runtime.  It is not
a replacement for a proof-assistant arithmetic kernel.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
from fractions import Fraction
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import mpmath as mp


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "FIN_P365_Independent_Check.json"


@dataclass(frozen=True)
class Interval:
    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @staticmethod
    def point(value: Fraction | int) -> "Interval":
        q = Fraction(value)
        return Interval(q, q)

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other: "Interval") -> "Interval":
        return self + (-other)

    def __mul__(self, other: "Interval") -> "Interval":
        values = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return Interval(min(values), max(values))

    def __truediv__(self, other: "Interval") -> "Interval":
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError
        reciprocal = Interval(
            min(1 / other.lo, 1 / other.hi),
            max(1 / other.lo, 1 / other.hi),
        )
        return self * reciprocal

    def scale(self, value: Fraction | int) -> "Interval":
        return self * Interval.point(value)

    def power(self, exponent: int) -> "Interval":
        result = Interval.point(1)
        for _ in range(exponent):
            result = result * self
        return result

    @property
    def mid(self) -> Fraction:
        return (self.lo + self.hi) / 2

    def strictly_inside(self, other: "Interval") -> bool:
        return other.lo < self.lo and self.hi < other.hi


def interval_sum(values: list[Interval]) -> Interval:
    result = Interval.point(0)
    for value in values:
        result = result + value
    return result


def inverse(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    n = len(matrix)
    augmented = [
        list(row) + [Fraction(int(i == j)) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = next(
            row for row in range(column, n) if augmented[row][column] != 0
        )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    left - factor * right
                    for left, right in zip(augmented[row], augmented[column])
                ]
    return [row[n:] for row in augmented]


def matrix_vector(
    matrix: list[list[Fraction]], vector: list[Interval]
) -> list[Interval]:
    return [
        interval_sum(
            [vector[column].scale(matrix[row][column]) for column in range(len(vector))]
        )
        for row in range(len(matrix))
    ]


def matrix_interval_product(
    left: list[list[Fraction]], right: list[list[Interval]]
) -> list[list[Interval]]:
    return [
        [
            interval_sum(
                [right[k][column].scale(left[row][k]) for k in range(len(right))]
            )
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def fifth_root_bracket(integer: int, digits: int = 60) -> tuple[Fraction, Fraction]:
    if integer == 0:
        return Fraction(0), Fraction(0)
    mp.mp.dps = digits + 40
    scale = 10**digits
    floor_value = int(mp.floor(mp.root(integer, 5) * scale))
    lower = Fraction(floor_value, scale)
    upper = Fraction(floor_value + 1, scale)
    assert lower**5 <= integer <= upper**5
    return lower, upper


def attenuation(order: int) -> Interval:
    if order == 0:
        return Interval.point(1)
    lower, upper = fifth_root_bracket(order)
    return Interval(
        1 / (Fraction(1) + upper**9),
        1 / (Fraction(1) + lower**9),
    )


def cos_interval(x: Fraction, terms: int = 42) -> Interval:
    value = Fraction(0)
    for n in range(terms + 1):
        term = x ** (2 * n) / math.factorial(2 * n)
        value += term if n % 2 == 0 else -term
    remainder = abs(x) ** (2 * terms + 2) / math.factorial(2 * terms + 2)
    return Interval(value - remainder, value + remainder)


def oscillatory(order: int) -> Interval:
    cosine = cos_interval(Fraction(743 * order, 4000) + Fraction(13, 80))
    if order == 0:
        return cosine
    lower, upper = fifth_root_bracket(order)
    return cosine / Interval(Fraction(1) + lower**9, Fraction(1) + upper**9)


def affine_substitution(
    coefficients: list[Fraction], left: Fraction, right: Fraction
) -> list[Fraction]:
    degree = len(coefficients) - 1
    width = right - left
    result = [Fraction(0)] * (degree + 1)
    for power, coefficient in enumerate(coefficients):
        for t_power in range(power + 1):
            result[t_power] += (
                coefficient
                * math.comb(power, t_power)
                * left ** (power - t_power)
                * width**t_power
            )
    return result


def to_bernstein(power: list[Fraction]) -> list[Fraction]:
    degree = len(power) - 1
    return [
        sum(
            (
                power[j] * Fraction(math.comb(i, j), math.comb(degree, j))
                for j in range(i + 1)
            ),
            Fraction(0),
        )
        for i in range(degree + 1)
    ]


def bernstein_range(
    coefficients: list[Fraction], subdivisions: int
) -> tuple[Fraction, Fraction]:
    lower: Fraction | None = None
    upper: Fraction | None = None
    for index in range(subdivisions):
        local = affine_substitution(
            coefficients,
            Fraction(index, subdivisions),
            Fraction(index + 1, subdivisions),
        )
        values = to_bernstein(local)
        lower = min(values) if lower is None else min(lower, *values)
        upper = max(values) if upper is None else max(upper, *values)
    assert lower is not None and upper is not None
    return lower, upper


def descriptor(value: Fraction) -> dict[str, Any]:
    numerator = value.numerator
    denominator = value.denominator
    numerator_bytes = abs(numerator).to_bytes(
        max(1, (abs(numerator).bit_length() + 7) // 8), "big"
    )
    denominator_bytes = denominator.to_bytes(
        max(1, (denominator.bit_length() + 7) // 8), "big"
    )
    digest = hashlib.sha256(
        (b"-" if numerator < 0 else b"+")
        + len(numerator_bytes).to_bytes(8, "big")
        + numerator_bytes
        + denominator_bytes
    ).hexdigest()
    return {
        "sha256": digest,
        "numerator_bits": abs(numerator).bit_length(),
        "denominator_bits": denominator.bit_length(),
    }


def parse_json_cell(cell: str) -> dict[str, Any]:
    return json.loads(cell)


def check_p351() -> dict[str, Any]:
    rows = list(
        csv.DictReader(
            (ROOT / "FIN_Programs_351_364_Envelope_Interval.csv").open(
                encoding="utf-8"
            )
        )
    )
    centers = [
        Fraction(str(value))
        for value in (
            0.23727834187306726,
            0.5481980169862998,
            0.8527659803199023,
            0.9418342033805291,
            0.3936855195357094,
            0.07118661177601973,
        )
    ]
    radius = Fraction(1, 10**12)
    box = [Interval(center - radius, center + radius) for center in centers]
    point = [Interval.point(center) for center in centers]
    moments = [attenuation(order) for order in range(7)]

    def system(variables: list[Interval]) -> list[Interval]:
        roots, weights = variables[:3], variables[3:]
        return [
            interval_sum(
                [weights[i] * roots[i].power(order) for i in range(3)]
            )
            - moments[order]
            for order in range(1, 7)
        ]

    def jacobian(variables: list[Interval]) -> list[list[Interval]]:
        roots, weights = variables[:3], variables[3:]
        return [
            [
                weights[i]
                * roots[i].power(order - 1)
                * Interval.point(order)
                for i in range(3)
            ]
            + [roots[i].power(order) for i in range(3)]
            for order in range(1, 7)
        ]

    c = inverse([[entry.mid for entry in row] for row in jacobian(point)])
    correction = matrix_vector(c, system(point))
    base = [Interval.point(centers[i]) - correction[i] for i in range(6)]
    cj = matrix_interval_product(c, jacobian(box))
    remainder_matrix = [
        [
            Interval.point(int(row == column)) - cj[row][column]
            for column in range(6)
        ]
        for row in range(6)
    ]
    delta = [Interval(-radius, radius) for _ in range(6)]
    remainder = [
        interval_sum(
            [remainder_matrix[row][column] * delta[column] for column in range(6)]
        )
        for row in range(6)
    ]
    image = [base[i] + remainder[i] for i in range(6)]
    inside = [image[i].strictly_inside(box[i]) for i in range(6)]
    identity_matches = []
    for row, value in zip(rows, image):
        identity_matches.extend(
            [
                descriptor(value.lo) == parse_json_cell(row["krawczyk_lower_exact"]),
                descriptor(value.hi) == parse_json_cell(row["krawczyk_upper_exact"]),
            ]
        )
    return {
        "all_components_strictly_inside": all(inside),
        "endpoint_identity_checks": len(identity_matches),
        "all_endpoint_identities_match": all(identity_matches),
    }


def check_p352() -> dict[str, Any]:
    archived = list(
        csv.DictReader(
            (ROOT / "FIN_Programs_351_364_Oscillatory_Interval.csv").open(
                encoding="utf-8"
            )
        )
    )
    source = list(
        csv.DictReader(
            (ROOT / "FIN_Programs_337_350_Oscillatory_Resource.csv").open(
                encoding="utf-8"
            )
        )
    )
    raw = [
        Fraction(row["value"])
        for row in source
        if row["row_type"] == "continuum_dual_polynomial"
    ]
    raw_lower, raw_upper = bernstein_range(raw, 4096)
    span = raw_upper - raw_lower
    safe = [value / span for value in raw]
    safe[0] -= raw_upper / span
    safe_lower, safe_upper = bernstein_range(safe, 4096)

    safe_rows = [
        row for row in archived if row["row_type"] == "safe_dual_coefficient"
    ]
    safe_matches = [
        descriptor(value) == parse_json_cell(row["exact_value"])
        for value, row in zip(safe, safe_rows)
    ]

    nodes = [
        Fraction(79, 4000),
        Fraction(80, 4000),
        Fraction(517, 4000),
        Fraction(518, 4000),
        Fraction(1180, 4000),
        Fraction(1181, 4000),
        Fraction(2107, 4000),
        Fraction(3256, 4000),
        Fraction(3257, 4000),
        Fraction(3767, 4000),
        Fraction(3768, 4000),
        Fraction(1),
    ]
    moments = [oscillatory(order) for order in range(12)]
    vandermonde = [[node**order for node in nodes] for order in range(12)]
    weights = matrix_vector(inverse(vandermonde), moments)
    weight_rows = [
        row for row in archived if row["row_type"] == "fixed_support_weight_interval"
    ]
    weight_matches = []
    signs = []
    for value, row in zip(weights, weight_rows):
        weight_matches.extend(
            [
                descriptor(value.lo) == parse_json_cell(row["weight_lower_exact"]),
                descriptor(value.hi) == parse_json_cell(row["weight_upper_exact"]),
            ]
        )
        signs.append(1 if value.lo > 0 else -1 if value.hi < 0 else 0)
    return {
        "safe_bernstein_range": [float(safe_lower), float(safe_upper)],
        "dual_is_globally_feasible": safe_lower >= -1 and safe_upper <= 0,
        "safe_coefficient_identity_checks": len(safe_matches),
        "all_safe_coefficient_identities_match": all(safe_matches),
        "weight_endpoint_identity_checks": len(weight_matches),
        "all_weight_endpoint_identities_match": all(weight_matches),
        "all_weight_signs_certified": 0 not in signs,
        "positive_weights": signs.count(1),
        "negative_weights": signs.count(-1),
    }


def main() -> None:
    result = {
        "checker": "P365 independent standard-library exact arithmetic",
        "imports_FIN_research_modules": False,
        "P351": check_p351(),
        "P352": check_p352(),
        "boundary": (
            "This is an independent implementation in the same Python "
            "runtime. It does not replace proof-assistant kernel checking."
        ),
    }
    result["all_checks_pass"] = bool(
        result["P351"]["all_components_strictly_inside"]
        and result["P351"]["all_endpoint_identities_match"]
        and result["P352"]["dual_is_globally_feasible"]
        and result["P352"]["all_safe_coefficient_identities_match"]
        and result["P352"]["all_weight_endpoint_identities_match"]
        and result["P352"]["all_weight_signs_certified"]
    )
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)
    if not result["all_checks_pass"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
