#!/usr/bin/env python3
"""Independent exact-rational enclosure engine for FIN ST82.

The proof layer uses only Python integers and Fraction arithmetic.  Binary64
linear algebra supplies centers and preconditioners, but every accepted
Krawczyk inequality is re-evaluated exactly after rationalization.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import root

from fin_st01_st15_research import N, strict_operator
from fin_st16_st27_research import final_memory_state


ROOT = Path(__file__).resolve().parent
CERTIFICATE = ROOT / "FIN_ST82_Independent_Rational_Certificate.json"


@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("invalid interval")

    @staticmethod
    def point(value) -> "I":
        item = value if isinstance(value, F) else F(value)
        return I(item, item)

    def __add__(self, other) -> "I":
        other = as_i(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self) -> "I":
        return I(-self.hi, -self.lo)

    def __sub__(self, other) -> "I":
        return self + (-as_i(other))

    def __rsub__(self, other) -> "I":
        return as_i(other) - self

    def __mul__(self, other) -> "I":
        other = as_i(other)
        values = [self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi]
        return I(min(values), max(values))

    __rmul__ = __mul__

    def reciprocal(self) -> "I":
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError("interval contains zero")
        return I(min(1 / self.lo, 1 / self.hi), max(1 / self.lo, 1 / self.hi))

    def __truediv__(self, other) -> "I":
        return self * as_i(other).reciprocal()

    def __rtruediv__(self, other) -> "I":
        return as_i(other) / self

    def square(self) -> "I":
        if self.lo <= 0 <= self.hi:
            return I(F(0), max(self.lo * self.lo, self.hi * self.hi))
        return I(min(self.lo * self.lo, self.hi * self.hi), max(self.lo * self.lo, self.hi * self.hi))

    def power(self, exponent: int) -> "I":
        if exponent < 0:
            return self.power(-exponent).reciprocal()
        if exponent == 0:
            return I.point(1)
        if exponent % 2 == 0:
            return self.square().power(exponent // 2)
        return self * self.power(exponent - 1)

    def midpoint(self) -> F:
        return (self.lo + self.hi) / 2

    def maxabs(self) -> F:
        return max(abs(self.lo), abs(self.hi))

    def width(self) -> F:
        return self.hi - self.lo

    def contains_interior(self, other: "I") -> bool:
        return self.lo < other.lo and other.hi < self.hi


def as_i(value) -> I:
    return value if isinstance(value, I) else I.point(value)


def float_fraction(value: float) -> F:
    return F(format(float(value), ".17g"))


def arctan_alternating(x: F, terms: int) -> I:
    total = F(0)
    for k in range(terms):
        total += (-1 if k % 2 else 1) * x ** (2 * k + 1) / (2 * k + 1)
    next_term = (-1 if terms % 2 else 1) * x ** (2 * terms + 1) / (2 * terms + 1)
    return I(min(total, total + next_term), max(total, total + next_term))


def pi_interval(terms: int = 90) -> I:
    a = arctan_alternating(F(1, 5), terms)
    b = arctan_alternating(F(1, 239), terms)
    return 16 * a - 4 * b


def trig_taylor(argument: I, cosine: bool, terms: int = 48) -> I:
    total = I.point(0)
    if cosine:
        for k in range(terms):
            total += ((-1) ** k) * argument.power(2 * k) / math.factorial(2 * k)
        order = 2 * terms
    else:
        for k in range(terms):
            total += ((-1) ** k) * argument.power(2 * k + 1) / math.factorial(2 * k + 1)
        order = 2 * terms + 1
    remainder = argument.maxabs() ** order / math.factorial(order)
    return I(total.lo - remainder, total.hi + remainder)


def reduce_pi_coefficient(value: F) -> F:
    # Reduce modulo 2 into [-1,1).
    quotient = (value + 1) // 2
    return value - 2 * quotient


def cos_pi(coefficient: F, pi: I) -> I:
    coefficient = reduce_pi_coefficient(coefficient)
    return trig_taylor(coefficient * pi, cosine=True)


def sin_pi(coefficient: F, pi: I) -> I:
    coefficient = reduce_pi_coefficient(coefficient)
    return trig_taylor(coefficient * pi, cosine=False)


def integer_nth_root_floor(value: int, n: int) -> int:
    low, high = 0, 1
    while high**n <= value:
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if middle**n <= value:
            low = middle
        else:
            high = middle
    return low


def fifth_root_interval(value: int, digits: int = 38) -> I:
    scale = 10**digits
    floor = integer_nth_root_floor(value * scale**5, 5)
    return I(F(floor, scale), F(floor + 1, scale))


def spectral_data() -> dict:
    pi = pi_interval()
    omega, phi = F(743, 4000), F(13, 80)
    weights: dict[int, I] = {}
    for distance in range(1, 7):
        numerator = trig_taylor(I.point(omega * distance + phi), cosine=True)
        denominator = I.point(1) + fifth_root_interval(distance).power(9)
        weights[distance] = numerator / denominator
    row_sum = 2 * sum((weights[d] for d in range(1, 6)), I.point(0)) + weights[6]

    lambdas = []
    for k in range(12):
        if k == 0:
            value = I.point(0)
        else:
            value = row_sum
            for distance in range(1, 6):
                value -= 2 * weights[distance] * cos_pi(F(k * distance, 6), pi)
            value -= weights[6] * ((-1) ** k)
        lambdas.append(value)

    poles, residues = [], []
    for k, multiplicity in zip([0, 1, 2], [1, 2, 2]):
        pole = row_sum + F(3, 20)
        pole -= 2 * weights[2] * cos_pi(F(k, 3), pi)
        pole -= 2 * weights[4] * cos_pi(F(2 * k, 3), pi)
        pole -= weights[6] * ((-1) ** k)
        poles.append(pole)
        real, imag = I.point(0), I.point(0)
        for j, distance in enumerate([1, 3, 5, 5, 3, 1]):
            entry = -weights[distance]
            real += entry * cos_pi(F(k * j, 3), pi)
            imag -= entry * sin_pi(F(k * j, 3), pi)
        residues.append(F(multiplicity, 6) * (real.square() + imag.square()))

    fzero = sum((residue / pole for pole, residue in zip(poles, residues)), I.point(0))
    memory, tvalues, svalues = [], [], []
    for lam in lambdas:
        memory.append(fzero - sum((residue / (lam + pole) for pole, residue in zip(poles, residues)), I.point(0)))
        tvalues.append(sum((residue / (lam + pole).power(2) for pole, residue in zip(poles, residues)), I.point(0)))
        svalues.append(sum((residue / (lam + pole).power(3) for pole, residue in zip(poles, residues)), I.point(0)))
    return {
        "pi": pi,
        "weights": weights,
        "row_sum": row_sum,
        "lambdas": lambdas,
        "poles": poles,
        "residues": residues,
        "memory": memory,
        "tvalues": tvalues,
        "svalues": svalues,
    }


def circulant_from_eigenvalues(values: list[I], pi: I) -> list[list[I]]:
    first = []
    for j in range(N):
        item = I.point(0)
        for k, value in enumerate(values):
            item += value * cos_pi(F(k * j, 6), pi) / N
        first.append(item)
    return [[first[(column - row) % N] for column in range(N)] for row in range(N)]


def point_matrix(values: np.ndarray) -> list[list[F]]:
    return [[float_fraction(values[row, column]) for column in range(values.shape[1])] for row in range(values.shape[0])]


def midpoint_matrix(values: list[list[I]]) -> np.ndarray:
    return np.asarray([[float(item.midpoint()) for item in row] for row in values], dtype=float)


def point_vector(values: np.ndarray) -> list[I]:
    return [I.point(float_fraction(value)) for value in values]


def imatvec(matrix: list[list[I]], vector: list[I]) -> list[I]:
    return [sum((matrix[row][column] * vector[column] for column in range(len(vector))), I.point(0)) for row in range(len(matrix))]


def pmat_imat(left: list[list[F]], right: list[list[I]]) -> list[list[I]]:
    return [
        [sum((left[row][k] * right[k][column] for k in range(len(right))), I.point(0)) for column in range(len(right[0]))]
        for row in range(len(left))
    ]


def pmat_ivec(matrix: list[list[F]], vector: list[I]) -> list[I]:
    return [sum((matrix[row][column] * vector[column] for column in range(len(vector))), I.point(0)) for row in range(len(matrix))]


def identity_minus(values: list[list[I]]) -> list[list[I]]:
    output = [[-item for item in row] for row in values]
    for index in range(len(output)):
        output[index][index] += 1
    return output


def centered_box(center: np.ndarray, radius: F) -> tuple[list[I], list[I]]:
    points = [float_fraction(value) for value in center]
    box = [I(value - radius, value + radius) for value in points]
    deltas = [I(-radius, radius) for _ in points]
    return box, deltas


def nonlinear_certificate(operator: list[list[I]], candidate: np.ndarray, radius: F = F(1, 10**12)) -> dict:
    midpoint = midpoint_matrix(operator)
    solution = root(
        lambda u: midpoint @ u + u - u**3,
        candidate,
        jac=lambda u: midpoint + np.eye(N) - np.diag(3 * u**2),
        options={"xtol": 1e-12, "maxfev": 5000},
    ).x
    center = point_vector(solution)
    box, delta = centered_box(solution, radius)
    jacobian_mid = midpoint + np.eye(N) - np.diag(3 * solution**2)
    inverse = point_matrix(np.linalg.inv(jacobian_mid))
    fcenter = imatvec(operator, center)
    for i in range(N):
        fcenter[i] += center[i] - center[i].power(3)
    correction = [-item for item in pmat_ivec(inverse, fcenter)]
    jacobian = [[item for item in row] for row in operator]
    for i in range(N):
        jacobian[i][i] += 1 - 3 * box[i].square()
    defect = identity_minus(pmat_imat(inverse, jacobian))
    image = [correction[i] + imatvec(defect, delta)[i] for i in range(N)]
    margins = [radius - item.maxabs() for item in image]
    if min(margins) <= 0:
        raise RuntimeError(f"nonlinear rational Krawczyk failed: {float(min(margins))}")
    return {
        "solution": solution,
        "box": box,
        "radius": radius,
        "minimum_margin": min(margins),
        "maximum_image_radius": max(item.maxabs() for item in image),
    }


def linear_certificate(matrix: list[list[I]], rhs: list[I], radius: F = F(1, 10**11)) -> dict:
    midpoint = midpoint_matrix(matrix)
    rhs_midpoint = np.asarray([float(item.midpoint()) for item in rhs])
    solution = np.linalg.solve(midpoint, rhs_midpoint)
    center = point_vector(solution)
    inverse = point_matrix(np.linalg.inv(midpoint))
    residual = [rhs[i] - imatvec(matrix, center)[i] for i in range(len(rhs))]
    correction = pmat_ivec(inverse, residual)
    defect = identity_minus(pmat_imat(inverse, matrix))
    delta = [I(-radius, radius) for _ in solution]
    defect_delta = imatvec(defect, delta)
    image = [correction[i] + defect_delta[i] for i in range(len(solution))]
    margins = [radius - item.maxabs() for item in image]
    if min(margins) <= 0:
        raise RuntimeError(f"linear rational Krawczyk failed: {float(min(margins))}")
    box = [I(float_fraction(value) - radius, float_fraction(value) + radius) for value in solution]
    return {
        "solution": solution,
        "box": box,
        "radius": radius,
        "minimum_margin": min(margins),
        "maximum_image_radius": max(item.maxabs() for item in image),
    }


def idot(left: list[I], right: list[I]) -> I:
    return sum((a * b for a, b in zip(left, right)), I.point(0))


def compact(interval: I) -> list[float]:
    return [float(interval.lo), float(interval.hi)]


def run_certificate() -> dict:
    data = spectral_data()
    operator_eigenvalues = [a + m for a, m in zip(data["lambdas"], data["memory"])]
    operator = circulant_from_eigenvalues(operator_eigenvalues, data["pi"])
    top = circulant_from_eigenvalues(data["tvalues"], data["pi"])
    second = circulant_from_eigenvalues(data["svalues"], data["pi"])
    _, strict_mid, _ = strict_operator()
    candidate = final_memory_state(strict_mid)[0]
    nonlinear = nonlinear_certificate(operator, candidate)
    state = nonlinear["box"]

    linearized = [[item for item in row] for row in operator]
    for i in range(N):
        linearized[i][i] += 1 - 3 * state[i].square()
    y = linear_certificate(linearized, state)
    v = imatvec(top, state)
    z = linear_certificate(linearized, v)
    su = imatvec(second, state)
    c0 = idot(state, y["box"])
    c1 = idot(state, z["box"])
    c2 = idot(v, z["box"])
    c3 = idot(state, su)
    coefficients = {
        "a": c2 + c3,
        "b": 2 * c1,
        "c": c0,
    }

    # Exact-rational monotone root enclosure.  A deliberately broad initial
    # bracket is bisected only when coefficient-interval signs are decisive.
    left, right = F("1.277"), F("1.279")

    def polynomial(point: F) -> I:
        return coefficients["a"] + coefficients["b"] * point + coefficients["c"] * point * point

    if polynomial(left).lo <= 0 or polynomial(right).hi >= 0:
        raise RuntimeError("independent collision root lacks initial sign bracket")
    for _ in range(50):
        middle = (left + right) / 2
        value = polynomial(middle)
        if value.lo > 0:
            left = middle
        elif value.hi < 0:
            right = middle
        else:
            break
    derivative = coefficients["b"] + 2 * coefficients["c"] * I(left, right)
    if derivative.hi >= 0:
        raise RuntimeError("independent collision polynomial is not decreasing")

    previous = json.loads((ROOT / "FIN_ST58_Full_Interval_Certificate.json").read_text(encoding="utf-8"))
    previous_interval = previous["collision_speed_interval"]
    overlap = not (float(right) < previous_interval[0] or float(left) > previous_interval[1])
    certificate = {
        "arithmetic": "exact Python Fraction intervals; Machin pi; alternating arctan; Taylor remainder; integer fifth-root bracketing; exact rational Krawczyk inequalities",
        "independent_of": "mpmath.iv and all ST58 serialized interval endpoints except final comparison",
        "pi_interval": compact(data["pi"]),
        "maximum_weight_interval_width": float(max(item.width() for item in data["weights"].values())),
        "strict_eigenvalue_intervals": [compact(item) for item in data["lambdas"]],
        "pole_intervals": [compact(item) for item in data["poles"]],
        "residue_intervals": [compact(item) for item in data["residues"]],
        "nonlinear_state_radius": float(nonlinear["radius"]),
        "nonlinear_minimum_inclusion_margin": float(nonlinear["minimum_margin"]),
        "linear_minimum_inclusion_margins": [float(y["minimum_margin"]), float(z["minimum_margin"])],
        "coefficient_intervals": {key: compact(value) for key, value in coefficients.items()},
        "collision_speed_interval": [float(left), float(right)],
        "collision_interval_width": float(right - left),
        "collision_derivative_interval": compact(derivative),
        "overlaps_ST58_interval": overlap,
        "proof_boundary": "The exact proof uses rationalized binary64 centers and inverse matrices only as candidate preconditioners; exact Fraction Krawczyk inclusion pays their approximation. It remains a source-code certificate, not a proof-assistant theorem or global Krein result.",
    }
    CERTIFICATE.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    return certificate


if __name__ == "__main__":
    result = run_certificate()
    print(json.dumps(result, indent=2, sort_keys=True))
