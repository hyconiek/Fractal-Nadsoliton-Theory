#!/usr/bin/env python3
"""FIN local research Programs P445--P447.

P445 derives an exact diagonal-support reduction of the two-slot reduced
phase/dephasing comb and closes the selected nonideal cell analytically.
P446 certifies the global optimum inside the palindromic three-use symmetric
code line and audits, without promoting, the full simplex frontier.
P447 propagates the complete P429 atom boxes through the detector-envelope
minimax Jordan sampler.

All inputs are local mathematical artifacts; no physical record is consumed.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import heapq
import itertools
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize, minimize_scalar

import fin_programs_411_427 as p411
import fin_programs_435_436_440 as p435


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_445_447"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
P445_PATH = ROOT / f"{PREFIX}_P445_Two_Slot_Reduction.csv"
P446_PATH = ROOT / f"{PREFIX}_P446_Palindromic_Certificate.csv"
P446_SIMPLEX_PATH = ROOT / f"{PREFIX}_P446_Full_Simplex_Audit.csv"
P447_PATH = ROOT / f"{PREFIX}_P447_Atom_Box_Propagation.csv"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p445_p447_exact_reduction_and_intervals.png"
P429_CSV = ROOT / "FIN_Programs_428_430_P429_Krawczyk.csv"

Interval = tuple[Fraction, Fraction]


def json_ready(value: Any) -> Any:
    if isinstance(value, Fraction):
        if max(value.numerator.bit_length(), value.denominator.bit_length()) <= 512:
            return str(value)
        return format(float(value), ".17g")
    if isinstance(value, np.generic):
        return value.item()
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
            writer.writerow({key: (
                json.dumps(json_ready(row.get(key)), ensure_ascii=False)
                if isinstance(row.get(key), (list, tuple, dict))
                else json_ready(row.get(key, ""))
            ) for key in keys})


def iv_add(left: Interval, right: Interval) -> Interval:
    return left[0] + right[0], left[1] + right[1]


def iv_neg(value: Interval) -> Interval:
    return -value[1], -value[0]


def iv_sub(left: Interval, right: Interval) -> Interval:
    return iv_add(left, iv_neg(right))


def iv_mul(left: Interval, right: Interval) -> Interval:
    values = [left[i] * right[j] for i in (0, 1) for j in (0, 1)]
    return min(values), max(values)


def iv_reciprocal(value: Interval) -> Interval:
    if value[0] <= 0 <= value[1]:
        raise ZeroDivisionError("interval contains zero")
    return min(1 / value[0], 1 / value[1]), max(1 / value[0], 1 / value[1])


def iv_div(left: Interval, right: Interval) -> Interval:
    return iv_mul(left, iv_reciprocal(right))


def iv_scale(value: Interval, scalar: Fraction | int) -> Interval:
    return iv_mul(value, (Fraction(scalar), Fraction(scalar)))


def iv_square(value: Interval) -> Interval:
    if value[0] <= 0 <= value[1]:
        return Fraction(0), max(value[0] ** 2, value[1] ** 2)
    return min(value[0] ** 2, value[1] ** 2), max(value[0] ** 2, value[1] ** 2)


def iv_abs(value: Interval) -> Interval:
    if value[0] <= 0 <= value[1]:
        return Fraction(0), max(abs(value[0]), abs(value[1]))
    return min(abs(value[0]), abs(value[1])), max(abs(value[0]), abs(value[1]))


def iv_sqrt(value: Interval, places: int = 24) -> Interval:
    if value[0] < 0:
        raise ValueError("negative interval passed to square root")
    return (
        p435.rational_sqrt_interval(value[0], places)[0],
        p435.rational_sqrt_interval(value[1], places)[1],
    )


def iv_sum(values: list[Interval]) -> Interval:
    result: Interval = (Fraction(0), Fraction(0))
    for value in values:
        result = iv_add(result, value)
    return result


# ---------------------------------------------------------------------------
# P445: exact two-slot diagonal-support comb reduction
# ---------------------------------------------------------------------------


def p445_reduced_distance(
    history_probabilities: np.ndarray,
    q: float,
    theta: float,
    normalizer_coherence: complex = 0j,
) -> float:
    """Closed two-slot tester distance after diagonal-support reduction."""

    if history_probabilities.shape != (4,) or np.any(history_probabilities < -1e-15):
        raise ValueError("four nonnegative history probabilities required")
    probabilities = np.maximum(history_probabilities, 0)
    probabilities /= probabilities.sum()
    a = 2 * q * abs(math.sin(theta))
    b = 2 * q**2 * abs(math.sin(2 * theta))
    t0, t1, t2, t3 = probabilities
    u = float(np.real(normalizer_coherence))
    squared = (
        a * a * (t0 + t3) * (t1 + t2)
        + b * b * t0 * t3
        - 4 * a * a * u * u
        + 2 * a * b * (t3 - t0) * u
    )
    return math.sqrt(max(0.0, squared))


def program_445() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    q = Fraction(4, 5)
    theta = "pi/8"
    sqrt2 = math.sqrt(2)
    a = 2 * float(q) * math.sin(math.pi / 8)
    b = 2 * float(q) ** 2 * math.sin(math.pi / 4)
    condition_margin = b * b - 2 * a * a

    # In exact radicals: 2 q^2 cos^2(theta)-1=(8 sqrt(2)-9)/25.
    exact_condition = "(8*sqrt(2)-9)/25 > 0, because 128 > 81"
    exact_half_distance = "8*sqrt(2)/25"
    optimal_half_distance = 8 * sqrt2 / 25
    optimal_success = (1 + optimal_half_distance) / 2

    rng = np.random.default_rng(445)
    random_defects = []
    for probabilities in rng.dirichlet(np.ones(4), size=128):
        radius = math.sqrt(min(
            probabilities[0] * probabilities[1],
            probabilities[2] * probabilities[3],
        )) * math.sqrt(rng.random())
        phase = rng.uniform(0, 2 * math.pi)
        coherence = radius * np.exp(1j * phase)
        # Direct four-dimensional support computation.
        j_plus = p435.qubit_phase_dephasing_choi(float(q), math.pi / 8, +1)
        j_minus = p435.qubit_phase_dephasing_choi(float(q), math.pi / 8, -1)
        delta = np.kron(j_plus, j_plus) - np.kron(j_minus, j_minus)
        support = np.array([0, 3, 12, 15])
        compressed = delta[np.ix_(support, support)]
        normalizer = np.zeros((4, 4), dtype=complex)
        normalizer[:2, :2] = [
            [probabilities[0], coherence],
            [np.conjugate(coherence), probabilities[1]],
        ]
        normalizer[2:, 2:] = [
            [probabilities[2], -coherence],
            [-np.conjugate(coherence), probabilities[3]],
        ]
        eigenvalues, eigenvectors = np.linalg.eigh(normalizer)
        root = (eigenvectors * np.sqrt(np.maximum(eigenvalues, 0))) @ eigenvectors.conj().T
        direct = 0.5 * np.sum(np.abs(np.linalg.eigvalsh(root @ compressed @ root)))
        random_defects.append(abs(direct - p445_reduced_distance(
            probabilities, float(q), math.pi / 8, coherence
        )))

    # The optimized variable is x=t00+t11.  For fixed x, t00*t11<=x^2/4.
    coefficient_a2 = a * a
    coefficient_b2_over4 = b * b / 4
    rows = []
    for x in np.linspace(0, 1, 21):
        value = math.sqrt(max(0.0, coefficient_a2 * x * (1 - x) + coefficient_b2_over4 * x * x))
        rows.append({
            "diagonal_history_mass_x": x,
            "optimized_half_distance_at_x": value,
            "endpoint_product_t00_t11": x * x / 4,
        })

    return ({
        "status": (
            "[Proven] exact diagonal-support characterization of every two-slot tester "
            "for the declared phase/dephasing comb and exact GHZ optimality at q=4/5, theta=pi/8"
        ),
        "history_simplex": (
            "t=(t00,t01,t10,t11), t>=0, sum(t)=1, with normalizer coherence "
            "c=u+iv satisfying |c|^2<=min(t00*t01,t10*t11)"
        ),
        "reduced_formula": (
            "D(t,c)^2=a^2(t00+t11)(t01+t10)+b^2*t00*t11-4a^2*u^2+"
            "2ab(t11-t00)u, a=2q sin(theta), b=2q^2 sin(2theta); Im(c) cancels"
        ),
        "square_completion": (
            "D^2=F(x)-4a^2(u-b(t11-t00)/(4a))^2, "
            "F(x)=a^2*x*(1-x)+(b^2/4)*x^2"
        ),
        "one_variable_upper": "D(t,c)^2<=F(x), x=t00+t11",
        "ghz_condition": "b^2 >= 2*a^2",
        "exact_condition_check": exact_condition,
        "condition_margin_float": condition_margin,
        "optimal_history_law": ["1/2", "0", "0", "1/2"],
        "optimal_half_distance_exact": exact_half_distance,
        "optimal_half_distance": optimal_half_distance,
        "optimal_success": optimal_success,
        "maximum_random_formula_defect": max(random_defects),
        "theorem_scope": (
            "All causally admissible two-slot testers for the declared reduced memoryless "
            "qubit channel; arbitrary intermediate memory is included by tester normalization."
        ),
        "boundary": (
            "The proof uses the diagonal Choi support and two slots. It does not extend "
            "automatically to three or more uses, noncommuting noise, unheralded loss, "
            "or the full twelve-mode FIN channel."
        ),
        "new_object": "O161 Diagonal-Support Comb Reduction",
    }, rows)


# ---------------------------------------------------------------------------
# P446: palindromic interval closure and full-simplex audit
# ---------------------------------------------------------------------------


INV_SQRT2 = iv_sqrt((Fraction(1, 2), Fraction(1, 2)), 28)
INV_SQRT3 = iv_sqrt((Fraction(1, 3), Fraction(1, 3)), 28)


def interval_basis(rows: int, columns: int, entries: dict[tuple[int, int], Interval]) -> list[list[Interval]]:
    result = [[(Fraction(0), Fraction(0)) for _ in range(columns)] for _ in range(rows)]
    for (row, column), value in entries.items():
        result[row][column] = value
    return result


U2 = interval_basis(4, 3, {
    (0, 0): (Fraction(1), Fraction(1)),
    (1, 1): INV_SQRT2, (2, 1): INV_SQRT2,
    (3, 2): (Fraction(1), Fraction(1)),
})
U3 = interval_basis(8, 4, {
    (0, 0): (Fraction(1), Fraction(1)),
    (1, 1): INV_SQRT3, (2, 1): INV_SQRT3, (4, 1): INV_SQRT3,
    (3, 2): INV_SQRT3, (5, 2): INV_SQRT3, (6, 2): INV_SQRT3,
    (7, 3): (Fraction(1), Fraction(1)),
})
V3 = interval_basis(8, 2, {
    (1, 0): INV_SQRT2, (2, 0): iv_neg(INV_SQRT2),
    (6, 1): INV_SQRT2, (5, 1): iv_neg(INV_SQRT2),
})

U2_NUM = np.zeros((4, 3))
U2_NUM[0, 0] = 1
U2_NUM[[1, 2], 1] = 1 / math.sqrt(2)
U2_NUM[3, 2] = 1
U3_NUM = np.zeros((8, 4))
U3_NUM[0, 0] = 1
U3_NUM[[1, 2, 4], 1] = 1 / math.sqrt(3)
U3_NUM[[3, 5, 6], 2] = 1 / math.sqrt(3)
U3_NUM[7, 3] = 1
V3_NUM = np.zeros((8, 2))
V3_NUM[1, 0] = 1 / math.sqrt(2)
V3_NUM[2, 0] = -1 / math.sqrt(2)
V3_NUM[6, 1] = 1 / math.sqrt(2)
V3_NUM[5, 1] = -1 / math.sqrt(2)


def interval_project(matrix: list[list[Interval]], basis: list[list[Interval]]) -> list[list[Interval]]:
    source_dimension = len(matrix)
    target_dimension = len(basis[0])
    result = [[(Fraction(0), Fraction(0)) for _ in range(target_dimension)] for _ in range(target_dimension)]
    for left in range(target_dimension):
        for right in range(target_dimension):
            terms = []
            for row in range(source_dimension):
                if basis[row][left] == (0, 0):
                    continue
                for column in range(source_dimension):
                    if basis[column][right] == (0, 0):
                        continue
                    terms.append(iv_mul(iv_mul(basis[row][left], matrix[row][column]), basis[column][right]))
            result[left][right] = iv_sum(terms)
    return result


def interval_skew_matrix(probabilities: list[Interval], survivors: int) -> list[list[Interval]]:
    uses = 3
    lost = uses - survivors
    dimension = 2**survivors
    coherence = Fraction(4, 5)
    result = [[(Fraction(0), Fraction(0)) for _ in range(dimension)] for _ in range(dimension)]
    for left in range(dimension):
        left_weight_kept = left.bit_count()
        for right in range(dimension):
            right_weight_kept = right.bit_count()
            hamming = (left ^ right).bit_count()
            difference = left_weight_kept - right_weight_kept
            amplitude_terms = []
            for lost_value in range(2**lost):
                lost_weight = lost_value.bit_count()
                left_weight = left_weight_kept + lost_weight
                right_weight = right_weight_kept + lost_weight
                radicand = iv_scale(
                    iv_mul(probabilities[left_weight], probabilities[right_weight]),
                    Fraction(1, math.comb(3, left_weight) * math.comb(3, right_weight)),
                )
                amplitude_terms.append(iv_sqrt(radicand, 24))
            result[left][right] = iv_scale(
                iv_mul(iv_sum(amplitude_terms), p435.p436_sine_interval(difference)),
                2 * coherence**hamming,
            )
    return result


def skew_three_distance(matrix: list[list[Interval]]) -> Interval:
    squared = iv_sum([iv_square(matrix[i][j]) for i in range(3) for j in range(i + 1, 3)])
    return iv_sqrt(squared, 22)


def skew_four_distance(matrix: list[list[Interval]]) -> Interval:
    squared = iv_sum([iv_square(matrix[i][j]) for i in range(4) for j in range(i + 1, 4)])
    pfaffian = iv_add(
        iv_sub(iv_mul(matrix[0][1], matrix[2][3]), iv_mul(matrix[0][2], matrix[1][3])),
        iv_mul(matrix[0][3], matrix[1][2]),
    )
    radicand = iv_add(squared, iv_scale(iv_abs(pfaffian), 2))
    return iv_sqrt(radicand, 22)


def p446_objective_interval(probabilities: list[Interval]) -> Interval:
    k1 = interval_skew_matrix(probabilities, 1)
    distance1 = iv_abs(k1[0][1])

    k2 = interval_skew_matrix(probabilities, 2)
    distance2 = skew_three_distance(interval_project(k2, U2))

    k3 = interval_skew_matrix(probabilities, 3)
    symmetric_distance = skew_four_distance(interval_project(k3, U3))
    standard_block = interval_project(k3, V3)
    distance3 = iv_add(symmetric_distance, iv_scale(iv_abs(standard_block[0][1]), 2))

    return iv_sum([
        iv_scale(distance1, Fraction(12, 125)),
        iv_scale(distance2, Fraction(48, 125)),
        iv_scale(distance3, Fraction(64, 125)),
    ])


def p446_numeric_objective(probabilities: np.ndarray) -> float:
    """Fast S3-block formula, independently checked against P418 matrices."""

    matrices = []
    for survivors in range(1, 4):
        dimension = 2**survivors
        lost = 3 - survivors
        matrix = np.zeros((dimension, dimension))
        for left in range(dimension):
            for right in range(dimension):
                wl0 = left.bit_count(); wr0 = right.bit_count()
                amplitude = 0.0
                for lost_value in range(2**lost):
                    lost_weight = lost_value.bit_count()
                    wl = wl0 + lost_weight; wr = wr0 + lost_weight
                    amplitude += math.sqrt(
                        max(0.0, probabilities[wl] * probabilities[wr])
                        / (math.comb(3, wl) * math.comb(3, wr))
                    )
                matrix[left, right] = (
                    2 * 0.8 ** ((left ^ right).bit_count())
                    * math.sin((2 * math.pi / 15) * (wl0 - wr0))
                    * amplitude
                )
        matrices.append(matrix)
    d1 = abs(matrices[0][0, 1])
    block2 = U2_NUM.T @ matrices[1] @ U2_NUM
    d2 = math.sqrt(sum(block2[i, j] ** 2 for i in range(3) for j in range(i + 1, 3)))
    block3 = U3_NUM.T @ matrices[2] @ U3_NUM
    entries = [block3[i, j] for i in range(4) for j in range(i + 1, 4)]
    pfaffian = block3[0, 1] * block3[2, 3] - block3[0, 2] * block3[1, 3] + block3[0, 3] * block3[1, 2]
    d3s = math.sqrt(sum(value * value for value in entries) + 2 * abs(pfaffian))
    standard = V3_NUM.T @ matrices[2] @ V3_NUM
    d3 = d3s + 2 * abs(standard[0, 1])
    return (12 * d1 + 48 * d2 + 64 * d3) / 125


# Directed binary64 interval arithmetic for the high-throughput one-dimensional
# cover.  Every primitive operation is expanded by one representable number.
FInterval = tuple[float, float]


def fdown(value: float) -> float:
    return float(np.nextafter(value, -np.inf))


def fup(value: float) -> float:
    return float(np.nextafter(value, np.inf))


def fiv_from_fraction(value: Interval) -> FInterval:
    low = float(value[0]); high = float(value[1])
    if Fraction.from_float(low) > value[0]:
        low = fdown(low)
    if Fraction.from_float(high) < value[1]:
        high = fup(high)
    return fdown(low), fup(high)


def fadd(left: FInterval, right: FInterval) -> FInterval:
    return fdown(left[0] + right[0]), fup(left[1] + right[1])


def fneg(value: FInterval) -> FInterval:
    return -value[1], -value[0]


def fsub(left: FInterval, right: FInterval) -> FInterval:
    return fadd(left, fneg(right))


def fmul(left: FInterval, right: FInterval) -> FInterval:
    values = [left[i] * right[j] for i in (0, 1) for j in (0, 1)]
    return fdown(min(values)), fup(max(values))


def fscale(value: FInterval, scalar: FInterval) -> FInterval:
    return fmul(value, scalar)


def fsquare(value: FInterval) -> FInterval:
    if value[0] <= 0 <= value[1]:
        return 0.0, fup(max(value[0] ** 2, value[1] ** 2))
    return fdown(min(value[0] ** 2, value[1] ** 2)), fup(max(value[0] ** 2, value[1] ** 2))


def fabsiv(value: FInterval) -> FInterval:
    if value[0] <= 0 <= value[1]:
        return 0.0, fup(max(abs(value[0]), abs(value[1])))
    return fdown(min(abs(value[0]), abs(value[1]))), fup(max(abs(value[0]), abs(value[1])))


def fsqrt(value: FInterval) -> FInterval:
    return fdown(math.sqrt(max(0.0, value[0]))), fup(math.sqrt(max(0.0, value[1])))


def fsum(values: list[FInterval]) -> FInterval:
    result = (0.0, 0.0)
    for value in values:
        result = fadd(result, value)
    return result


Q_FI = fiv_from_fraction((Fraction(4, 5), Fraction(4, 5)))
SINES_FI = {difference: fiv_from_fraction(p435.p436_sine_interval(difference)) for difference in range(-3, 4)}
INV_SQRT2_FI = fiv_from_fraction(INV_SQRT2)
INV_SQRT3_FI = fiv_from_fraction(INV_SQRT3)


def float_basis(rows: int, columns: int, entries: dict[tuple[int, int], FInterval]) -> list[list[FInterval]]:
    result = [[(0.0, 0.0) for _ in range(columns)] for _ in range(rows)]
    for key, value in entries.items():
        result[key[0]][key[1]] = value
    return result


U2_FI = float_basis(4, 3, {(0, 0): (1, 1), (1, 1): INV_SQRT2_FI, (2, 1): INV_SQRT2_FI, (3, 2): (1, 1)})
U3_FI = float_basis(8, 4, {(0, 0): (1, 1), (1, 1): INV_SQRT3_FI, (2, 1): INV_SQRT3_FI, (4, 1): INV_SQRT3_FI, (3, 2): INV_SQRT3_FI, (5, 2): INV_SQRT3_FI, (6, 2): INV_SQRT3_FI, (7, 3): (1, 1)})
V3_FI = float_basis(8, 2, {(1, 0): INV_SQRT2_FI, (2, 0): fneg(INV_SQRT2_FI), (6, 1): INV_SQRT2_FI, (5, 1): fneg(INV_SQRT2_FI)})


def fproject(matrix: list[list[FInterval]], basis: list[list[FInterval]]) -> list[list[FInterval]]:
    source = len(matrix); target = len(basis[0])
    result = [[(0.0, 0.0) for _ in range(target)] for _ in range(target)]
    for a in range(target):
        for b in range(target):
            terms = []
            for i in range(source):
                if basis[i][a] == (0.0, 0.0):
                    continue
                for j in range(source):
                    if basis[j][b] != (0.0, 0.0):
                        terms.append(fmul(fmul(basis[i][a], matrix[i][j]), basis[j][b]))
            result[a][b] = fsum(terms)
    return result


def fskew(probabilities: list[FInterval], survivors: int) -> list[list[FInterval]]:
    dimension = 2**survivors; lost = 3-survivors
    result = [[(0.0, 0.0) for _ in range(dimension)] for _ in range(dimension)]
    for left in range(dimension):
        for right in range(dimension):
            wl0=left.bit_count(); wr0=right.bit_count(); h=(left^right).bit_count()
            terms=[]
            for lv in range(2**lost):
                z=lv.bit_count(); wl=wl0+z; wr=wr0+z
                rad=fscale(fmul(probabilities[wl],probabilities[wr]),(1/(math.comb(3,wl)*math.comb(3,wr)),1/(math.comb(3,wl)*math.comb(3,wr))))
                terms.append(fsqrt(rad))
            qh=(Q_FI[0]**h,Q_FI[1]**h)
            result[left][right]=fscale(fmul(fsum(terms),SINES_FI[wl0-wr0]),fscale(qh,(2,2)))
    return result


def fobjective(probabilities: list[FInterval]) -> FInterval:
    m1=fskew(probabilities,1);d1=fabsiv(m1[0][1])
    b2=fproject(fskew(probabilities,2),U2_FI)
    d2=fsqrt(fsum([fsquare(b2[i][j]) for i in range(3) for j in range(i+1,3)]))
    b3=fproject(fskew(probabilities,3),U3_FI)
    sq=fsum([fsquare(b3[i][j]) for i in range(4) for j in range(i+1,4)])
    pf=fadd(fsub(fmul(b3[0][1],b3[2][3]),fmul(b3[0][2],b3[1][3])),fmul(b3[0][3],b3[1][2]))
    d3s=fsqrt(fadd(sq,fscale(fabsiv(pf),(2,2))))
    std=fproject(fskew(probabilities,3),V3_FI)
    d3=fadd(d3s,fscale(fabsiv(std[0][1]),(2,2)))
    return fsum([fscale(d1,(12/125,12/125)),fscale(d2,(48/125,48/125)),fscale(d3,(64/125,64/125))])


def palindromic_probabilities(interval_a: Interval) -> list[Interval]:
    interval_b = (Fraction(1, 2) - interval_a[1], Fraction(1, 2) - interval_a[0])
    return [interval_a, interval_b, interval_b, interval_a]


def palindromic_branch_certificate(tolerance: float = 1e-3) -> dict[str, Any]:
    def probabilities(interval_a: FInterval) -> list[FInterval]:
        interval_b = (fdown(0.5 - interval_a[1]), fup(0.5 - interval_a[0]))
        return [interval_a, interval_b, interval_b, interval_a]

    initial: FInterval = (0.0, 0.5)
    initial_value = fobjective(probabilities(initial))
    heap: list[tuple[float, int, FInterval, FInterval]] = []
    counter = 0
    heapq.heappush(heap, (-float(initial_value[1]), counter, initial, initial_value))
    # Import the already certified P436 rational feasible value as incumbent.
    incumbent_lower = 0.46327828319203235
    incumbent_a = 56923 / 250000
    terminal: list[tuple[FInterval, FInterval]] = []
    evaluated = 0
    while heap:
        negative_upper, _, interval_a, value_interval = heapq.heappop(heap)
        if value_interval[1] <= incumbent_lower + tolerance:
            terminal.append((interval_a, value_interval))
            continue
        midpoint = (interval_a[0] + interval_a[1]) / 2
        midpoint_value = fobjective(probabilities((midpoint, midpoint)))
        evaluated += 1
        if midpoint_value[0] > incumbent_lower:
            incumbent_lower = midpoint_value[0]
            incumbent_a = midpoint
        if interval_a[1] - interval_a[0] <= 1e-12:
            terminal.append((interval_a, value_interval))
            continue
        for child in ((interval_a[0], midpoint), (midpoint, interval_a[1])):
            child_value = fobjective(probabilities(child))
            counter += 1
            if child_value[1] > incumbent_lower:
                heapq.heappush(heap, (-float(child_value[1]), counter, child, child_value))
            else:
                terminal.append((child, child_value))
        if evaluated > 200000:
            raise RuntimeError("palindromic interval certificate exceeded its finite budget")

    global_upper = max([value[1] for _, value in terminal] + [incumbent_lower])
    surviving = [
        (box, value) for box, value in terminal
        if value[1] >= incumbent_lower
    ]
    maximizer_hull = (
        min(box[0] for box, _ in surviving),
        max(box[1] for box, _ in surviving),
    )
    return {
        "incumbent_a": incumbent_a,
        "incumbent_lower": incumbent_lower,
        "global_upper": global_upper,
        "optimality_gap": global_upper - incumbent_lower,
        "maximizer_hull": maximizer_hull,
        "evaluated_boxes": evaluated,
        "terminal_boxes": len(terminal),
    }


def full_simplex_scout() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    uses = 3
    theta = 2 * math.pi / 15
    coherence = survival = 0.8
    rng = np.random.default_rng(446)
    starts = [
        np.full(4, 0.25),
        np.array([1 / 8, 3 / 8, 3 / 8, 1 / 8]),
        np.array([0.5, 0, 0, 0.5]),
    ]
    starts.extend(rng.dirichlet(np.ones(4)) for _ in range(21))
    best_value = -1.0
    best_vector = starts[0]
    rows = []
    for index, start in enumerate(starts):
        fit = minimize(
            lambda vector: -p446_numeric_objective(vector),
            start,
            method="SLSQP",
            bounds=[(0, 1)] * 4,
            constraints={"type": "eq", "fun": lambda vector: vector.sum() - 1},
            options={"ftol": 1e-13, "maxiter": 500},
        )
        vector = np.maximum(fit.x, 0)
        vector /= vector.sum()
        value = p446_numeric_objective(vector)
        rows.append({
            "start": index,
            "value": value,
            "probabilities": vector.tolist(),
            "reversal_defect": float(np.max(np.abs(vector - vector[::-1]))),
            "solver_success": bool(fit.success),
        })
        if value > best_value:
            best_value, best_vector = value, vector

    # Cheap finite falsification atlas, not an interval cover.
    hessian_maxima = []

    def reduced_value(first_three: np.ndarray) -> float:
        vector = np.r_[first_three, 1 - first_three.sum()]
        return p446_numeric_objective(vector)

    for vector in rng.dirichlet(np.full(4, 2.0), size=16):
        point = vector[:3]
        step = min(1e-5, float(vector.min()) / 10)
        hessian = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                ei = np.zeros(3); ej = np.zeros(3)
                ei[i] = step; ej[j] = step
                hessian[i, j] = (
                    reduced_value(point + ei + ej)
                    - reduced_value(point + ei - ej)
                    - reduced_value(point - ei + ej)
                    + reduced_value(point - ei - ej)
                ) / (4 * step * step)
        hessian_maxima.append(float(np.linalg.eigvalsh((hessian + hessian.T) / 2)[-1]))

    return ({
        "status": (
            "[Strong evidence] all deterministic multistarts converge to the palindromic "
            "certificate neighborhood and sampled interior Hessians are negative; "
            "[Open] rigorous global full-simplex concavity or interval cover"
        ),
        "starts": len(starts),
        "best_value": best_value,
        "best_vector": best_vector.tolist(),
        "maximum_best_reversal_defect": float(np.max(np.abs(best_vector - best_vector[::-1]))),
        "maximum_sampled_hessian_eigenvalue": max(hessian_maxima),
        "sampled_hessian_points": len(hessian_maxima),
        "boundary": (
            "SLSQP convergence and sampled negative Hessians are falsification evidence, not "
            "a global concavity theorem. The palindromic interval proof cannot be promoted "
            "to the full three-dimensional simplex."
        ),
    }, rows)


def program_446() -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    certificate = palindromic_branch_certificate()
    scout, scout_rows = full_simplex_scout()
    rows = [{
        "incumbent_a": certificate["incumbent_a"],
        "incumbent_lower": certificate["incumbent_lower"],
        "global_upper": certificate["global_upper"],
        "optimality_gap": certificate["optimality_gap"],
        "maximizer_lower": certificate["maximizer_hull"][0],
        "maximizer_upper": certificate["maximizer_hull"][1],
        "evaluated_boxes": certificate["evaluated_boxes"],
        "terminal_boxes": certificate["terminal_boxes"],
    }]
    result = {
        "status": (
            "[Computer-assisted proof] a global 1e-3 palindromic-line upper/lower certificate "
            "from directed interval representation reduction; [Strong evidence / Open] full simplex"
        ),
        "palindromic_family": "p(a)=(a,1/2-a,1/2-a,a), 0<=a<=1/2",
        "palindromic_certificate": certificate,
        "full_simplex_audit": scout,
        "proof_method": (
            "S3 irreducible block reduction (4D symmetric plus two repeated 2D standard "
            "blocks), outward-rounded binary64 intervals, and exhaustive one-dimensional branch-and-bound"
        ),
        "boundary": (
            "Globality is proved only on the palindromic line. The full symmetric Hamming-sector "
            "simplex remains open despite multistart and sampled-Hessian evidence."
        ),
        "new_object": "O162 Palindromic Heralded-Code Global Upper Certificate",
    }
    return result, rows, scout_rows


# ---------------------------------------------------------------------------
# P447: full P429 box propagation through O160
# ---------------------------------------------------------------------------


def p429_atom_boxes() -> tuple[list[Interval], list[Interval]]:
    values: dict[str, Interval] = {}
    with P429_CSV.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["row_type"] != "certified_variable":
                continue
            name = row["variable"]
            if name.startswith("x") or name.startswith("w"):
                values[name] = (Fraction(row["box_lower"]), Fraction(row["box_upper"]))
    nodes = [values[f"x{index}"] for index in range(6)] + [(Fraction(1), Fraction(1))]
    weights = [values[f"w{index}"] for index in range(7)]
    return nodes, weights


def iv_power_nonnegative(value: Interval, exponent: int) -> Interval:
    return value[0] ** exponent, value[1] ** exponent


def program_447() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    nodes, weights = p429_atom_boxes()
    radii = []
    hbars = []
    absolute_weights = []
    for node, weight in zip(nodes, weights):
        radius_squared = iv_sum([iv_power_nonnegative(node, 2 * order) for order in range(12)])
        radii.append(iv_sqrt(radius_squared, 30))
        efficiency_low = iv_add((Fraction(13, 20), Fraction(13, 20)), iv_scale(node, Fraction(1, 4)))
        dark_upper = iv_sub((Fraction(3, 100), Fraction(3, 100)), iv_scale(node, Fraction(1, 50)))
        hbar = iv_add(
            iv_reciprocal(efficiency_low),
            iv_div(iv_mul(dark_upper, iv_sub((Fraction(1), Fraction(1)), dark_upper)), iv_square(efficiency_low)),
        )
        hbars.append(hbar)
        absolute_weights.append(iv_abs(weight))

    scores_baseline = absolute_weights
    scores_p422 = [iv_mul(weight, radius) for weight, radius in zip(absolute_weights, radii)]
    scores_p447 = [
        iv_mul(iv_mul(weight, radius), iv_sqrt(hbar, 30))
        for weight, radius, hbar in zip(absolute_weights, radii, hbars)
    ]

    def normalized_intervals(scores: list[Interval]) -> list[Interval]:
        denominator = iv_sum(scores)
        return [iv_div(score, denominator) for score in scores]

    q_baseline = normalized_intervals(scores_baseline)
    q_p422 = normalized_intervals(scores_p422)
    q_p447 = normalized_intervals(scores_p447)

    # Interval moments and worst MSE coefficients.
    moments = []
    for order in range(12):
        moments.append(iv_sum([
            iv_mul(weight, iv_power_nonnegative(node, order))
            for weight, node in zip(weights, nodes)
        ]))
    moment_norm_squared = iv_sum([iv_square(moment) for moment in moments])

    def mse_interval(probabilities: list[Interval]) -> Interval:
        second_terms = []
        for weight, radius, hbar, probability in zip(weights, radii, hbars, probabilities):
            numerator = iv_mul(iv_mul(iv_square(weight), iv_square(radius)), hbar)
            second_terms.append(iv_div(numerator, probability))
        return iv_sub(iv_sum(second_terms), moment_norm_squared)

    mse_baseline = mse_interval(q_baseline)
    mse_p422 = mse_interval(q_p422)
    mse_p447 = mse_interval(q_p447)
    reduction_baseline = iv_sub((Fraction(1), Fraction(1)), iv_div(mse_p447, mse_baseline))
    reduction_p422 = iv_sub((Fraction(1), Fraction(1)), iv_div(mse_p447, mse_p422))

    rows = []
    for index in range(7):
        rows.append({
            "atom": index,
            "node_lower": nodes[index][0], "node_upper": nodes[index][1],
            "weight_lower": weights[index][0], "weight_upper": weights[index][1],
            "P447_probability_lower": q_p447[index][0],
            "P447_probability_upper": q_p447[index][1],
            "P447_probability_width": q_p447[index][1] - q_p447[index][0],
            "P422_probability_lower": q_p422[index][0],
            "P422_probability_upper": q_p422[index][1],
            "worst_second_moment_factor_lower": hbars[index][0],
            "worst_second_moment_factor_upper": hbars[index][1],
        })

    return ({
        "status": (
            "[Computer-assisted proof] the complete P429 isolating boxes propagate through "
            "the declared O160 detector envelope without changing positivity, normalization, "
            "or the ordering of the three worst-MSE coefficients"
        ),
        "P447_probability_intervals": q_p447,
        "maximum_probability_width": max(value[1] - value[0] for value in q_p447),
        "mse_absolute_weight_interval": mse_baseline,
        "mse_P422_interval": mse_p422,
        "mse_P447_interval": mse_p447,
        "reduction_vs_absolute_weight_interval": reduction_baseline,
        "reduction_vs_P422_interval": reduction_p422,
        "strict_order_certified": (
            mse_p447[1] < mse_p422[0] and mse_p422[1] < mse_baseline[0]
        ),
        "proof_method": (
            "exact P429 rational boxes, rational interval arithmetic, and integer square-root enclosures"
        ),
        "boundary": (
            "This propagates mathematical atom uncertainty only. The detector nuisance envelope "
            "remains supplied and synthetic; no calibration or event record is generated."
        ),
        "new_object": "O163 Certified Detector-Allocation Tube",
    }, rows)


def make_figure(p445_rows: list[dict[str, Any]], p446: dict[str, Any], p447_rows: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.2))
    axes[0].plot(
        [row["diagonal_history_mass_x"] for row in p445_rows],
        [row["optimized_half_distance_at_x"] for row in p445_rows],
        "o-", color="#2563eb",
    )
    axes[0].axvline(1, color="#dc2626", linestyle="--", label="GHZ optimum")
    axes[0].set_title("P445 exact two-slot reduction")
    axes[0].set_xlabel("x = t00 + t11")
    axes[0].set_ylabel("half distance")
    axes[0].legend(fontsize=8)

    certificate = p446["palindromic_certificate"]
    lo, hi = [float(Fraction(value)) for value in certificate["maximizer_hull"]]
    grid = np.linspace(max(0, lo - 0.035), min(0.5, hi + 0.035), 120)
    values = [p446_numeric_objective(np.array([a, 0.5-a, 0.5-a, a])) for a in grid]
    axes[1].plot(grid, values, color="#059669")
    axes[1].axvspan(lo, hi, color="#f59e0b", alpha=0.35, label="certified maximizer hull")
    axes[1].set_title("P446 palindromic global certificate")
    axes[1].set_xlabel("a in (a,1/2-a,1/2-a,a)")
    axes[1].set_ylabel("heralded half distance")
    axes[1].legend(fontsize=8)

    atoms = [row["atom"] for row in p447_rows]
    lower = np.array([float(Fraction(row["P447_probability_lower"])) for row in p447_rows])
    upper = np.array([float(Fraction(row["P447_probability_upper"])) for row in p447_rows])
    centers = (lower + upper) / 2
    axes[2].errorbar(atoms, centers, yerr=(upper-lower)/2, fmt="o-", capsize=3, color="#7c3aed")
    axes[2].set_title("P447 certified allocation tube")
    axes[2].set_xlabel("P429 atom")
    axes[2].set_ylabel("minimax probability")
    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)


def main() -> None:
    p445_result, p445_rows = program_445()
    p446_result, p446_rows, p446_scout_rows = program_446()
    p447_result, p447_rows = program_447()
    payload = {
        "metadata": {
            "programs": ["P445", "P446", "P447"],
            "date": "2026-08-01",
            "local_only": True,
            "external_physical_evidence": False,
            "selector_discharged": False,
            "dimensional_source_exported": False,
            "legacy_strict_bridge_complete": False,
        },
        "P445": p445_result,
        "P446": p446_result,
        "P447": p447_result,
    }
    RESULTS_PATH.write_text(json.dumps(json_ready(payload), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_csv(P445_PATH, p445_rows)
    write_csv(P446_PATH, p446_rows)
    write_csv(P446_SIMPLEX_PATH, p446_scout_rows)
    write_csv(P447_PATH, p447_rows)
    write_csv(SUMMARY_PATH, [
        {"program": "P445", "status": p445_result["status"], "new_object": p445_result["new_object"], "boundary": p445_result["boundary"]},
        {"program": "P446", "status": p446_result["status"], "new_object": p446_result["new_object"], "boundary": p446_result["boundary"]},
        {"program": "P447", "status": p447_result["status"], "new_object": p447_result["new_object"], "boundary": p447_result["boundary"]},
    ])
    make_figure(p445_rows, p446_result, p447_rows)
    print(json.dumps({
        "p445_two_slot_half_distance": p445_result["optimal_half_distance"],
        "p446_palindromic_gap": float(Fraction(p446_result["palindromic_certificate"]["optimality_gap"])),
        "p446_full_simplex_best": p446_result["full_simplex_audit"]["best_value"],
        "p447_max_probability_width": float(Fraction(p447_result["maximum_probability_width"])),
        "p447_strict_order": p447_result["strict_order_certified"],
    }, indent=2))


if __name__ == "__main__":
    main()
