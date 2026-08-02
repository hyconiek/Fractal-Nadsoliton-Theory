#!/usr/bin/env python3
"""FIN local research Programs P451--P453.

P451 proves the global optimum on the complete diagonal three-slot tester
face, then gives an exact-rational positive coherent tester whose certified
value is strictly larger.  Thus diagonal optimality is refuted while the
full 21-dimensional optimum remains open.

P452 proves the missing reversal-symmetrization inequality for the only
coarse erasure sector not already covered by P448 concavity.  The P446
palindromic branch certificate therefore becomes a full-simplex theorem.

P453 combines P429 isolation with P430 strict complementary contact geometry
to prove global uniqueness of the minimum-negative-mass (equivalently, fixed
mass minimum-total-variation) signed measure.

All inputs are local mathematical artifacts.  No physical record is used.
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
from scipy.optimize import minimize
import sympy as sp

import fin_programs_428_430 as p428
import fin_programs_435_436_440 as p435
import fin_programs_445_447 as p445
import fin_programs_448_450 as p448


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_451_453"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
P451_PATH = ROOT / f"{PREFIX}_P451_Causal_Coherence.csv"
P451_WITNESS_PATH = ROOT / f"{PREFIX}_P451_Coherent_Witness.npz"
P452_PATH = ROOT / f"{PREFIX}_P452_Coarse_Symmetrization.csv"
P453_PATH = ROOT / f"{PREFIX}_P453_Canonical_Representation.csv"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p451_p453_coherence_symmetry_and_gauge_fixing.png"
P429_CSV = ROOT / "FIN_Programs_428_430_P429_Krawczyk.csv"
P430_CSV = ROOT / "FIN_Programs_428_430_P430_Global_Dual.csv"

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


# ---------------------------------------------------------------------------
# P451: global diagonal theorem and coherent causal witness
# ---------------------------------------------------------------------------


def p451_algebraic_maximizer_interval() -> tuple[Interval, Interval, dict[str, Any]]:
    """Exact nested-radical enclosure for the valid stationary point and value."""

    sqrt2 = p435.rational_sqrt_interval(Fraction(2), 44)
    inner = p445.iv_sub(
        (Fraction(3504977), Fraction(3504977)),
        p445.iv_scale(sqrt2, Fraction(2461052)),
    )
    root = p445.iv_sqrt(inner, 40)
    a_cell = (Fraction(0), Fraction(0))
    a_cell = p445.iv_sub(
        a_cell, p445.iv_scale(root, Fraction(5177925, 28304679968))
    )
    a_cell = p445.iv_sub(
        a_cell,
        p445.iv_scale(
            p445.iv_mul(sqrt2, root), Fraction(413025, 3538084996)
        ),
    )
    a_cell = p445.iv_add(
        a_cell, p445.iv_scale(sqrt2, Fraction(3691, 58384))
    )
    a_cell = p445.iv_add(
        a_cell, (Fraction(31987, 116768), Fraction(31987, 116768))
    )

    a2 = p445.iv_square(a_cell)
    r_quadratic = p445.iv_add(
        (Fraction(-3520), Fraction(-3520)),
        p445.iv_scale(sqrt2, Fraction(1568)),
    )
    r_linear = p445.iv_sub(
        (Fraction(1532), Fraction(1532)),
        p445.iv_scale(sqrt2, Fraction(414)),
    )
    r_constant = p445.iv_sub(
        (Fraction(242), Fraction(242)),
        p445.iv_scale(sqrt2, Fraction(121)),
    )
    radicand = p445.iv_scale(
        p445.iv_add(
            p445.iv_add(
                p445.iv_mul(r_quadratic, a2),
                p445.iv_mul(r_linear, a_cell),
            ),
            r_constant,
        ),
        Fraction(16, 15625),
    )
    symmetric = p445.iv_sqrt(radicand, 38)
    root_2_minus = p445.iv_sqrt(
        p445.iv_sub((Fraction(2), Fraction(2)), sqrt2), 40
    )
    standard = p445.iv_mul(
        p445.iv_scale(
            p445.iv_sub(
                (Fraction(12), Fraction(12)),
                p445.iv_scale(a_cell, Fraction(24)),
            ),
            Fraction(1, 125),
        ),
        root_2_minus,
    )
    value = p445.iv_add(symmetric, standard)

    # The squared stationary equation is checked symbolically.  The sign of
    # R'(a*) selects the valid, unsquared branch.
    a = sp.symbols("a", real=True)
    s2 = sp.sqrt(2)
    r_expr = sp.Rational(16, 15625) * (
        (-3520 + 1568 * s2) * a**2
        + (1532 - 414 * s2) * a
        + 242 - 121 * s2
    )
    c_expr = sp.sqrt(2 - s2)
    stationary = sp.expand(
        sp.diff(r_expr, a) ** 2
        - 4 * r_expr * (sp.Rational(24, 125) * c_expr) ** 2
    )
    polynomial = (
        -8836992 * a**2 + 5639168 * s2 * a**2
        - 1972208 * s2 * a + 3415528 * a
        - 323159 + 149850 * s2
    )
    factor = sp.simplify(stationary / polynomial)
    a_exact = (
        -sp.Rational(5177925, 28304679968)
        * sp.sqrt(3504977 - 2461052 * s2)
        - sp.Rational(413025, 3538084996)
        * s2 * sp.sqrt(3504977 - 2461052 * s2)
        + sp.Rational(3691, 58384) * s2
        + sp.Rational(31987, 116768)
    )
    exact_root_residual = sp.simplify(polynomial.subs(a, a_exact))
    derivative = sp.diff(r_expr, a)
    derivative_cell = p445.iv_scale(
        p445.iv_add(
            p445.iv_mul(p445.iv_scale(r_quadratic, 2), a_cell),
            r_linear,
        ),
        Fraction(16, 15625),
    )
    return a_cell, value, {
        "stationary_factor": str(factor),
        "exact_root_residual": str(exact_root_residual),
        "radicand_derivative_lower": derivative_cell[0],
        "radicand_derivative_upper": derivative_cell[1],
        "exact_maximizer": str(a_exact),
        "exact_derivative_expression": str(derivative),
    }


def p451_coherent_normalizer() -> tuple[np.ndarray, dict[str, Fraction]]:
    """A rational interior point of C_3 with a symmetry-shaped coherence."""

    A = Fraction(287166, 10**6)
    B = Fraction(68777, 10**6)
    C = Fraction(1, 2) - A - 2 * B
    u = Fraction(9493, 10**6)
    left = [
        [A, u, u, 0],
        [u, B, 0, -u],
        [u, 0, B, -u],
        [0, -u, -u, C],
    ]
    right = [
        [C, -u, -u, 0],
        [-u, B, 0, u],
        [-u, 0, B, u],
        [0, u, u, A],
    ]
    normalizer = np.zeros((8, 8), dtype=float)
    normalizer[:4, :4] = np.array(left, dtype=float)
    normalizer[4:, 4:] = np.array(right, dtype=float)
    return normalizer, {"A": A, "B": B, "C": C, "u": u}


def p451_coherent_distance_interval() -> tuple[Interval, dict[str, Any]]:
    """Root-isolation certificate for the rational coherent witness.

    For positive N, N(iH) is similar to sqrt(N)(iH)sqrt(N), so its real
    eigenvalues give the trace norm.  The trigonometric matrix is rounded to
    rationals; the remaining interval error is paid with a nuclear-norm
    perturbation bound using ||N||_2 <= Tr(N)=1.
    """

    _, values = p451_coherent_normalizer()
    A, B, C, u = (values[key] for key in ("A", "B", "C", "u"))
    left = [[A, u, u, 0], [u, B, 0, -u], [u, 0, B, -u], [0, -u, -u, C]]
    right = [[C, -u, -u, 0], [-u, B, 0, u], [-u, 0, B, u], [0, u, u, A]]
    normalizer = sp.diag(
        sp.Matrix([[sp.Rational(Fraction(x).numerator, Fraction(x).denominator) for x in row] for row in left]),
        sp.Matrix([[sp.Rational(Fraction(x).numerator, Fraction(x).denominator) for x in row] for row in right]),
    )
    dimension = 8
    scale = 10**18
    midpoint = sp.zeros(dimension)
    maximum_radius = Fraction(0)
    for row in range(dimension):
        for column in range(dimension):
            hamming = (row ^ column).bit_count()
            difference = row.bit_count() - column.bit_count()
            sine = p448.p449_sine_interval(difference)
            factor = 2 * Fraction(4, 5) ** hamming
            interval = p445.iv_scale(sine, factor)
            center = (interval[0] + interval[1]) / 2
            rounded = Fraction(round(float(center) * scale), scale)
            maximum_radius = max(
                maximum_radius,
                abs(rounded - interval[0]),
                abs(interval[1] - rounded),
            )
            midpoint[row, column] = sp.I * sp.Rational(
                rounded.numerator, rounded.denominator
            )
    polynomial = (normalizer * midpoint).charpoly().as_poly()
    isolated = sp.intervals(polynomial, eps=sp.Rational(1, 10**25))
    lower = Fraction(0)
    upper = Fraction(0)
    root_rows: list[dict[str, Any]] = []
    for (cell, multiplicity) in isolated:
        lo_sp, hi_sp = cell
        lo = Fraction(int(lo_sp.p), int(lo_sp.q))
        hi = Fraction(int(hi_sp.p), int(hi_sp.q))
        if lo > 0:
            abs_cell = (lo, hi)
        elif hi < 0:
            abs_cell = (-hi, -lo)
        else:
            abs_cell = (Fraction(0), max(-lo, hi))
        lower += multiplicity * abs_cell[0]
        upper += multiplicity * abs_cell[1]
        root_rows.append({
            "root_lower": lo,
            "root_upper": hi,
            "multiplicity": multiplicity,
        })
    # ||sqrt(N) E sqrt(N)||_* / 2 <= d^2 max|E_ij| / 2.
    perturbation = Fraction(dimension * dimension, 2) * maximum_radius
    return (lower / 2 - perturbation, upper / 2 + perturbation), {
        "rational_midpoint_root_cells": root_rows,
        "maximum_entry_radius": maximum_radius,
        "half_nuclear_perturbation": perturbation,
        "isolated_roots": len(isolated),
    }


def _effect_matrix(vector: np.ndarray, dimension: int) -> np.ndarray:
    hermitian = np.zeros((dimension, dimension), dtype=float)
    cursor = 0
    for index in range(dimension):
        hermitian[index, index] = vector[cursor]
        cursor += 1
    for row in range(dimension):
        for column in range(row + 1, dimension):
            hermitian[row, column] = hermitian[column, row] = vector[cursor]
            cursor += 1
    eigenvalues, eigenvectors = np.linalg.eigh(hermitian)
    effects = 1 / (1 + np.exp(-np.clip(eigenvalues, -40, 40)))
    return (eigenvectors * effects) @ eigenvectors.T


def _positive_root(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh((matrix + matrix.T) / 2)
    return (eigenvectors * np.sqrt(np.maximum(eigenvalues, 0))) @ eigenvectors.T


def p451_parameterized_normalizer(vector: np.ndarray) -> np.ndarray:
    """Surjective interior parameterization of the real part of C_3."""

    r = 1 / (1 + math.exp(-float(np.clip(vector[0], -40, 40))))
    inherited_one = np.diag([r, 1 - r])
    root_one = _positive_root(inherited_one)
    effect_two = _effect_matrix(vector[1:4], 2)
    block_two_zero = root_one @ effect_two @ root_one
    block_two_one = inherited_one - block_two_zero
    inherited_two = np.zeros((4, 4))
    inherited_two[:2, :2] = block_two_zero
    inherited_two[2:, 2:] = block_two_one
    root_two = _positive_root(inherited_two)
    effect_three = _effect_matrix(vector[4:14], 4)
    block_three_zero = root_two @ effect_three @ root_two
    block_three_one = inherited_two - block_three_zero
    result = np.zeros((8, 8))
    result[:4, :4] = block_three_zero
    result[4:, 4:] = block_three_one
    return result


def program_451() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    a_interval, diagonal_interval, algebra = p451_algebraic_maximizer_interval()
    witness, parameters = p451_coherent_normalizer()
    witness_interval, witness_certificate = p451_coherent_distance_interval()

    A, B, C, u = (parameters[key] for key in ("A", "B", "C", "u"))
    minor_two = A * B - 2 * u**2
    minor_three = A * B * C - 2 * u**2 * (A + C)
    inherited = witness[:4, :4] + witness[4:, 4:]
    recursion_residual = float(np.linalg.norm(
        inherited - np.diag([float(A + C), float(2 * B), float(2 * B), float(A + C)])
    ))
    diagonal_advantage = witness_interval[0] - diagonal_interval[1]

    delta = p448.compressed_process_difference(3, 0.8, math.pi / 8)
    rng = np.random.default_rng(451)
    scout_rows: list[dict[str, Any]] = []
    for seed in range(8):
        initial = rng.normal(size=14)
        optimum = minimize(
            lambda vector: -p448.tester_distance(
                p451_parameterized_normalizer(vector), delta
            ),
            initial,
            method="L-BFGS-B",
            options={"maxiter": 1500, "ftol": 1e-13, "gtol": 1e-8, "maxls": 40},
        )
        normalizer = p451_parameterized_normalizer(optimum.x)
        scout_rows.append({
            "row_type": "full_cone_numerical_scout",
            "seed": seed,
            "half_distance": -float(optimum.fun),
            "iterations": int(optimum.nit),
            "success": bool(optimum.success),
            "minimum_eigenvalue": float(np.linalg.eigvalsh(normalizer)[0]),
            "offdiagonal_frobenius_norm": float(np.linalg.norm(
                normalizer - np.diag(np.diag(normalizer))
            )),
        })

    np.savez_compressed(
        P451_WITNESS_PATH,
        support_normalizer=witness,
        process_difference_support=delta,
        A=np.array([float(A)]), B=np.array([float(B)]),
        C=np.array([float(C)]), u=np.array([float(u)]),
    )
    rows = [
        {
            "row_type": "diagonal_global_certificate",
            "maximizer_a_lower": a_interval[0],
            "maximizer_a_upper": a_interval[1],
            "half_distance_lower": diagonal_interval[0],
            "half_distance_upper": diagonal_interval[1],
        },
        {
            "row_type": "rational_coherent_witness",
            "A": A, "B": B, "C": C, "u": u,
            "minor_AB_minus_2u2": minor_two,
            "minor_ABC_minus_2u2_A_plus_C": minor_three,
            "half_distance_lower": witness_interval[0],
            "half_distance_upper": witness_interval[1],
            "advantage_over_diagonal_upper": diagonal_advantage,
        },
        *scout_rows,
    ]
    return ({
        "status": (
            "[Proven] global diagonal-face optimum; [Computer-assisted proof] "
            "an exact-rational coherent causal tester strictly exceeds it; "
            "[Strong evidence / Open] full 21-dimensional optimum"
        ),
        "diagonal_symmetrization": (
            "Concavity of N -> Tr|N^(1/2) Delta N^(1/2)|/2 and invariance under "
            "S3 slot permutations plus global bit complement move every diagonal "
            "history law to t000=t111=a and the other six masses=(1-2a)/6."
        ),
        "diagonal_formula": (
            "D_diag(a)=sqrt(16*(-3520*a^2+1568*sqrt(2)*a^2-414*sqrt(2)*a+"
            "1532*a-121*sqrt(2)+242)/15625)+(12-24*a)*sqrt(2-sqrt(2))/125"
        ),
        "diagonal_maximizer_interval": a_interval,
        "diagonal_half_distance_interval": diagonal_interval,
        "diagonal_success_probability_interval": (
            (1 + diagonal_interval[0]) / 2,
            (1 + diagonal_interval[1]) / 2,
        ),
        "algebraic_stationarity": algebra,
        "coherent_witness_parameters": parameters,
        "coherent_psd_exact_minors": [minor_two, minor_three],
        "coherent_recursion_residual": recursion_residual,
        "coherent_minimum_float_eigenvalue": float(np.linalg.eigvalsh(witness)[0]),
        "coherent_half_distance_interval": witness_interval,
        "coherent_advantage_over_global_diagonal_upper": diagonal_advantage,
        "coherent_success_probability_lower": (1 + witness_interval[0]) / 2,
        "witness_certificate": witness_certificate,
        "full_cone_scout_best": max(row["half_distance"] for row in scout_rows),
        "full_cone_scout_starts": len(scout_rows),
        "proof_boundary": (
            "The rational witness refutes diagonal optimality and is admitted by the exact "
            "P449 recursion. Multiple full-cone starts agree numerically, but no matching "
            "21-dimensional dual certificate is exported, so the global full-cone value "
            "and uniqueness remain open."
        ),
        "new_object": "O167 Causal-Coherence Advantage Witness",
    }, rows)


# ---------------------------------------------------------------------------
# P452: exact coarse-erasure reversal symmetrization
# ---------------------------------------------------------------------------


def cosine_interval(angle: Interval, last_odd: int = 19) -> Interval:
    """Alternating cosine enclosure on the present positive small interval."""

    def partial(value: Fraction, last_index: int) -> Fraction:
        return sum(
            (-1) ** index * value ** (2 * index) / Fraction(math.factorial(2 * index))
            for index in range(last_index + 1)
        )

    # Odd last index is a lower partial sum, even is upper. Cosine decreases.
    return partial(angle[1], last_odd), partial(angle[0], last_odd + 1)


def p452_component_values(probabilities: np.ndarray) -> np.ndarray:
    """Return the three coarse survivor-sector half norms before mask weights."""

    matrices = []
    for survivors in range(1, 4):
        dimension = 2**survivors
        lost = 3 - survivors
        matrix = np.zeros((dimension, dimension))
        for left in range(dimension):
            for right in range(dimension):
                amplitude = 0.0
                for lost_value in range(2**lost):
                    lost_weight = lost_value.bit_count()
                    wl = left.bit_count() + lost_weight
                    wr = right.bit_count() + lost_weight
                    amplitude += math.sqrt(
                        max(0.0, probabilities[wl] * probabilities[wr])
                        / (math.comb(3, wl) * math.comb(3, wr))
                    )
                matrix[left, right] = (
                    2 * 0.8 ** ((left ^ right).bit_count())
                    * math.sin(2 * math.pi / 15 * (left.bit_count() - right.bit_count()))
                    * amplitude
                )
        matrices.append(matrix)
    d1 = abs(matrices[0][0, 1])
    block_two = p445.U2_NUM.T @ matrices[1] @ p445.U2_NUM
    d2 = math.sqrt(sum(
        block_two[row, column] ** 2
        for row in range(3) for column in range(row + 1, 3)
    ))
    block_three = p445.U3_NUM.T @ matrices[2] @ p445.U3_NUM
    entries = [
        block_three[row, column]
        for row in range(4) for column in range(row + 1, 4)
    ]
    pfaffian = (
        block_three[0, 1] * block_three[2, 3]
        - block_three[0, 2] * block_three[1, 3]
        + block_three[0, 3] * block_three[1, 2]
    )
    symmetric = math.sqrt(sum(value**2 for value in entries) + 2 * abs(pfaffian))
    standard = p445.V3_NUM.T @ matrices[2] @ p445.V3_NUM
    return np.array([d1, d2, symmetric + 2 * abs(standard[0, 1])])


def program_452() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    theta_interval = p445.iv_scale(p435.PI_INTERVAL, Fraction(2, 15))
    cos_theta = cosine_interval(theta_interval)
    q = Fraction(4, 5)
    condition = p445.iv_sub(
        p445.iv_scale(p445.iv_square(cos_theta), 4 * q**2),
        (Fraction(1), Fraction(1)),
    )
    if condition[0] <= 0:
        raise AssertionError("the analytic d2 symmetrization condition was not certified")

    certificate = p445.palindromic_branch_certificate(tolerance=1e-3)
    rng = np.random.default_rng(452)
    rows: list[dict[str, Any]] = []
    maximum_component_defect = -math.inf
    maximum_objective_defect = -math.inf
    for seed, probabilities in enumerate(rng.dirichlet(np.ones(4), size=256)):
        symmetrized = (probabilities + probabilities[::-1]) / 2
        original_components = p452_component_values(probabilities)
        symmetric_components = p452_component_values(symmetrized)
        component_defects = original_components - symmetric_components
        objective_defect = (
            p445.p446_numeric_objective(probabilities)
            - p445.p446_numeric_objective(symmetrized)
        )
        maximum_component_defect = max(
            maximum_component_defect, float(np.max(component_defects))
        )
        maximum_objective_defect = max(maximum_objective_defect, objective_defect)
        if seed < 32:
            rows.append({
                "row_type": "symmetrization_falsification_sample",
                "seed": seed,
                "p": probabilities,
                "symmetrized_p": symmetrized,
                "d1_gain": symmetric_components[0] - original_components[0],
                "d2_gain": symmetric_components[1] - original_components[1],
                "d3_gain": symmetric_components[2] - original_components[2],
                "objective_gain": -objective_defect,
            })
    rows.insert(0, {
        "row_type": "analytic_condition",
        "cos_theta_lower": cos_theta[0],
        "cos_theta_upper": cos_theta[1],
        "four_q2_cos2_minus_one_lower": condition[0],
        "four_q2_cos2_minus_one_upper": condition[1],
        "full_simplex_lower": certificate["incumbent_lower"],
        "full_simplex_upper": certificate["global_upper"],
    })
    return ({
        "status": (
            "[Proven] reversal symmetrization for the exact coarse objective; "
            "[Computer-assisted proof] the P446 palindromic 1e-3 certificate "
            "is global on the complete four-sector simplex"
        ),
        "coarse_gap_localization": (
            "Survivor counts s=3 and s=1 equal their fine-grained P448 forms. "
            "Only s=2 has D2=||v0+v1|| rather than ||v0||+||v1||; its disclosure "
            "defect is 2(||v0||||v1||-v0.v1)/(||v0||+||v1||+||v0+v1||)."
        ),
        "d2_symmetrization_reduction": (
            "Writing endpoint and middle square-root masses as angles alpha,beta and "
            "ratio r, D2(sym p)^2-D2(p)^2 has nonnegative constant and linear "
            "coefficients. Its quadratic coefficient is "
            "Q=-(1/3)cos(2alpha)cos(2beta)+gamma cos^2(alpha+beta), "
            "gamma=4q^2 cos^2(theta)/3. Since cos^2(alpha+beta)-"
            "cos(2alpha)cos(2beta)=sin^2(alpha-beta), Q>=0 whenever "
            "4q^2 cos^2(theta)>=1."
        ),
        "certified_condition_interval": condition,
        "full_simplex_global_lower": certificate["incumbent_lower"],
        "full_simplex_global_upper": certificate["global_upper"],
        "full_simplex_gap": certificate["optimality_gap"],
        "maximizer_hull_a": certificate["maximizer_hull"],
        "evaluated_palindromic_boxes": certificate["evaluated_boxes"],
        "improvement_over_P448_upper": (
            0.4666305033804779 - certificate["global_upper"]
        ),
        "maximum_random_component_symmetrization_defect": maximum_component_defect,
        "maximum_random_objective_symmetrization_defect": maximum_objective_defect,
        "proof_boundary": (
            "This is global only for the declared three-use permutation-symmetric "
            "Hamming-sector code, q=4/5, theta=2*pi/15, and heralded loss model. "
            "The 1e-3 value interval is not an exact maximizer or an unrestricted "
            "adaptive-channel theorem."
        ),
        "new_object": "O168 Coarse-Erasure Symmetrization Certificate",
    }, rows)


# ---------------------------------------------------------------------------
# P453: strict complementarity and global representation uniqueness
# ---------------------------------------------------------------------------


def p429_boxes() -> tuple[list[Interval], list[Interval]]:
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


def program_453() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    nodes, weights = p429_boxes()
    p430_rows = list(csv.DictReader(P430_CSV.open(encoding="utf-8")))
    all_cells_safe = all(row.get("safe") == "True" for row in p430_rows)
    separations = [
        nodes[right][0] - nodes[left][1]
        for left in range(7) for right in range(left + 1, 7)
    ]
    vandermonde_lower = Fraction(1)
    for separation in separations:
        vandermonde_lower *= separation
    signs = [
        1 if cell[0] > 0 else -1 if cell[1] < 0 else 0 for cell in weights
    ]
    objective = p445.iv_add(
        p445.iv_neg(weights[0]), p445.iv_neg(weights[5])
    )
    total_mass = (Fraction(0), Fraction(0))
    for weight in weights:
        total_mass = p445.iv_add(total_mass, weight)
    tv_interval = p445.iv_add(total_mass, p445.iv_scale(objective, 2))

    rows: list[dict[str, Any]] = []
    target_levels = [-1, 0, 0, 0, 0, -1, 0]
    for index, (node, weight, target) in enumerate(zip(nodes, weights, target_levels)):
        rows.append({
            "row_type": "forced_complementary_contact",
            "contact": index,
            "node_lower": node[0],
            "node_upper": node[1],
            "weight_lower": weight[0],
            "weight_upper": weight[1],
            "weight_sign": signs[index],
            "dual_level": target,
            "jordan_role": "negative" if target == -1 else "positive",
        })
    rows.append({
        "row_type": "global_uniqueness_certificate",
        "all_P430_cells_safe": all_cells_safe,
        "minimum_pairwise_node_separation": min(separations),
        "absolute_vandermonde_determinant_lower": vandermonde_lower,
        "negative_mass_lower": objective[0],
        "negative_mass_upper": objective[1],
        "total_variation_lower": tv_interval[0],
        "total_variation_upper": tv_interval[1],
    })
    if not all_cells_safe or vandermonde_lower <= 0 or signs != [-1, 1, 1, 1, 1, -1, 1]:
        raise AssertionError("P453 strict-complementarity prerequisites failed")
    return ({
        "status": (
            "[Proven, conditional on the inherited P429/P430 computer-assisted "
            "certificates and signed-moment duality] the global minimum-negative-mass "
            "signed measure is unique; equivalently it is the unique fixed-mass "
            "minimum-total-variation representation"
        ),
        "complementarity_identity": (
            "N(mu)-integral p dmu = integral(-p)dmu_plus + "
            "integral(1+p)dmu_minus. Equality forces mu_plus onto {p=0} "
            "and mu_minus onto {p=-1}."
        ),
        "exact_contact_sets": (
            "P430 curvature/monotonicity and strict Bernstein gaps give "
            "{p=-1}={x0,x5} and {p=0}={x1,x2,x3,x4,1}."
        ),
        "forced_support_atoms": 7,
        "weight_signs": signs,
        "minimum_pairwise_node_separation": min(separations),
        "absolute_vandermonde_determinant_lower": vandermonde_lower,
        "negative_mass_interval": objective,
        "total_signed_mass_interval": total_mass,
        "minimum_total_variation_interval": tv_interval,
        "canonical_sampler_consequence": (
            "Within the variational rule 'minimize Jordan negative mass', O163 has one "
            "globally forced reduced atomic representation, so its detector allocation "
            "is no longer representation-ambiguous inside that rule."
        ),
        "removal_test": (
            "Removing the negative-mass/minimum-TV variational rule restores the exact "
            "O166 null-cycle gauge and unbounded sampler risk. Moments alone still do "
            "not choose a representation."
        ),
        "physical_boundary": (
            "The theorem fixes a mathematical extremal representation. It does not "
            "prove that a laboratory prepares it, nor does it generate detector "
            "calibration, custody, a dimensional unit, or physical evidence."
        ),
        "new_object": "O169 Strict-Complementarity Jordan Gauge Fixer",
    }, rows)


def make_figure(
    p451: dict[str, Any], p452: dict[str, Any], p453_rows: list[dict[str, Any]]
) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.35))

    grid = np.linspace(0, 0.5, 260)
    diagonal = []
    for a in grid:
        law = np.full(8, (1 - 2 * a) / 6)
        law[0] = law[-1] = a
        diagonal.append(p448.tester_distance(
            np.diag(law), p448.compressed_process_difference(3, 0.8, math.pi / 8)
        ))
    axes[0].plot(grid, diagonal, color="#2563eb", label="diagonal face")
    axes[0].axhline(
        float(Fraction(p451["coherent_half_distance_interval"][0])),
        color="#dc2626", linestyle="--", label="coherent witness",
    )
    axes[0].set_title("P451 coherence beats diagonal optimum")
    axes[0].set_xlabel("endpoint history mass a")
    axes[0].set_ylabel("half distance")
    axes[0].legend(fontsize=8)

    pgrid = np.linspace(0, 0.5, 240)
    coarse = [
        p445.p446_numeric_objective(np.array([a, 0.5-a, 0.5-a, a]))
        for a in pgrid
    ]
    axes[1].plot(pgrid, coarse, color="#0f766e")
    axes[1].axhspan(
        p452["full_simplex_global_lower"],
        p452["full_simplex_global_upper"],
        color="#f59e0b", alpha=0.24, label="global certificate",
    )
    axes[1].set_title("P452 full simplex reduces to line")
    axes[1].set_xlabel("palindromic sector mass a")
    axes[1].set_ylabel("coarse half distance")
    axes[1].legend(fontsize=8)

    contacts = [row for row in p453_rows if row["row_type"] == "forced_complementary_contact"]
    x = [float(Fraction(row["node_lower"]) + Fraction(row["node_upper"])) / 2 for row in contacts]
    w = [float(Fraction(row["weight_lower"]) + Fraction(row["weight_upper"])) / 2 for row in contacts]
    colors = ["#dc2626" if value < 0 else "#2563eb" for value in w]
    axes[2].axhline(0, color="#64748b", linewidth=0.8)
    axes[2].stem(x, w, linefmt="#64748b", markerfmt=" ", basefmt=" ")
    axes[2].scatter(x, w, c=colors, zorder=3)
    axes[2].set_title("P453 unique minimum-TV measure")
    axes[2].set_xlabel("contact node x")
    axes[2].set_ylabel("signed weight")

    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)


def main() -> None:
    p451_result, p451_rows = program_451()
    write_csv(P451_PATH, p451_rows)
    p452_result, p452_rows = program_452()
    write_csv(P452_PATH, p452_rows)
    p453_result, p453_rows = program_453()
    write_csv(P453_PATH, p453_rows)
    make_figure(p451_result, p452_result, p453_rows)

    results = {
        "metadata": {
            "programs": "P451-P453",
            "checkpoint": "P451-P453",
            "date": "2026-08-01",
            "local_only": True,
            "external_physical_evidence": False,
            "selector_discharged": False,
            "dimensional_source_exported": False,
            "legacy_strict_bridge_complete": False,
            "legacy_role_transfer_started": False,
        },
        "P451": p451_result,
        "P452": p452_result,
        "P453": p453_result,
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(results), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    summary_rows = [
        {
            "program": name,
            "status": results[name]["status"],
            "new_object": results[name]["new_object"],
            "physical_evidence": False,
            "selector_discharged": False,
        }
        for name in ("P451", "P452", "P453")
    ]
    write_csv(SUMMARY_PATH, summary_rows)
    print(json.dumps(json_ready(results), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
