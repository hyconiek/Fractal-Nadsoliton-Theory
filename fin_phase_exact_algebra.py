#!/usr/bin/env python3
"""Shared exact algebra for FIN Programs P485--P487.

This module performs only deterministic symbolic construction and exact
rational interval evaluation.  It does not launch an unbounded Groebner job.
"""

from __future__ import annotations

from fractions import Fraction
import itertools
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
CERTIFICATE_PATH = ROOT / "FIN_Program_480_Standalone_Certificate.json"
Interval = tuple[Fraction, Fraction]


def polynomial_system() -> dict[str, Any]:
    """Construct the portable exact 13-variable Riccati/KKT system."""

    A, B, u = sp.symbols("A B u", real=True)
    L, a, b, c, d, e, f, g, h, i = sp.symbols(
        "L a b c d e f g h i", real=True
    )
    s1, s2, s3 = sp.symbols("s1 s2 s3", real=True)
    variables = (A, B, u, L, a, b, c, d, e, f, g, h, i)
    C0 = sp.Rational(1, 2)-A-2*B
    left = sp.Matrix([[A,u,u,0],[u,B,0,-u],[u,0,B,-u],[0,-u,-u,C0]])
    right = sp.Matrix([[C0,-u,-u,0],[-u,B,0,u],[-u,0,B,u],[0,u,u,A]])
    normalizer = sp.diag(left, right)
    x2 = sp.Matrix([[L,a,a,c],[a,L,b,a],[a,b,L,a],[c,a,a,L]])
    cross = sp.Matrix([[d,e,e,f],[g,h,h,e],[g,h,h,e],[i,g,g,d]])
    x3 = x2.row_join(cross).col_join(cross.T.row_join(x2))
    q = sp.Rational(4, 5)
    sine = {0:0, 1:s1, 2:s2, 3:s3, -1:-s1, -2:-s2, -3:-s3}
    delta = sp.zeros(8)
    for row in range(8):
        for column in range(8):
            exponent = (row ^ column).bit_count()
            difference = row.bit_count()-column.bit_count()
            delta[row, column] = 2*sp.I*q**exponent*sine[difference]
    residual = sp.expand(x3*normalizer*x3-delta*normalizer*delta/4)
    groups: list[tuple[sp.Expr, list[tuple[int, int]]]] = []
    for row in range(8):
        for column in range(row, 8):
            expression = sp.expand(residual[row, column])
            for index, (representative, positions) in enumerate(groups):
                if sp.expand(expression-representative) == 0:
                    positions.append((row, column))
                    groups[index] = (representative, positions)
                    break
            else:
                groups.append((expression, [(row, column)]))
    if len(groups) != 13:
        raise AssertionError(f"expected 13 residual orbits, found {len(groups)}")
    equations = tuple(group[0] for group in groups)
    return {
        "variables": variables, "sines": (s1,s2,s3),
        "normalizer": normalizer, "x3": x3, "delta": delta,
        "equations": equations, "groups": [positions for _, positions in groups],
    }


def six_dimensional_symmetry_basis() -> sp.Matrix:
    root_two = sp.sqrt(2)
    basis = sp.zeros(8, 6)
    for row, value in ((0,1/root_two),(7,1/root_two)): basis[row,0] = value
    for row in (1,2,5,6): basis[row,1] = sp.Rational(1,2)
    for row, value in ((3,1/root_two),(4,1/root_two)): basis[row,2] = value
    for row, value in ((0,1/root_two),(7,-1/root_two)): basis[row,3] = value
    for row, value in ((1,sp.Rational(1,2)),(2,sp.Rational(1,2)),
                       (5,-sp.Rational(1,2)),(6,-sp.Rational(1,2))):
        basis[row,4] = value
    for row, value in ((3,1/root_two),(4,-1/root_two)): basis[row,5] = value
    return basis


def exact_representation_audit(system: dict[str, Any]) -> dict[str, Any]:
    transform = six_dimensional_symmetry_basis()
    reduced_n = sp.simplify(transform.T*system["normalizer"]*transform)
    reduced_x = sp.simplify(transform.T*system["x3"]*transform)
    reduced_k = sp.simplify(transform.T*(system["delta"]/sp.I)*transform)
    return {
        "N_plus": reduced_n[:3,:3], "X_plus": reduced_x[:3,:3],
        "X_minus": reduced_x[3:,3:], "C": reduced_k[:3,3:],
    }


def json_ready(value: Any) -> Any:
    if isinstance(value, sp.MatrixBase):
        return json_ready(value.tolist())
    if isinstance(value, sp.Basic):
        return str(value)
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def exact_context() -> dict[str, Any]:
    """Build the exact one-generator parity/Riccati context."""

    system = polynomial_system()
    representation = exact_representation_audit(system)
    alpha = sp.symbols("alpha", real=True)
    s1, s2, s3 = system["sines"]
    sine_substitution = {
        s1: alpha,
        s2: 1 - 2 * alpha**2,
        s3: 3 * alpha - 4 * alpha**3,
        sp.sqrt(2): 2 - 4 * alpha**2,
    }
    specialize = lambda value: sp.expand(value.subs(sine_substitution))
    n0 = specialize(representation["N_plus"])
    x_plus = specialize(representation["X_plus"])
    x_minus = specialize(representation["X_minus"])
    cross = specialize(representation["C"])
    plus_residual = sp.expand(x_plus*n0*x_plus-cross*n0*cross.T/4)
    minus_residual = sp.expand(x_minus*n0*x_minus-cross.T*n0*cross/4)
    active_equations = [
        sp.expand(matrix[row, column])
        for matrix in (plus_residual, minus_residual)
        for row in range(3)
        for column in range(row, 3)
    ]
    if len(active_equations) != 12:
        raise AssertionError("P485/P487 require twelve active parity equations")

    A, B, u, L, a, b, c, d, e, f, g, h, i = system["variables"]
    branch_equation = sp.expand(125*(L-b)-36*alpha)
    alpha_polynomial = 8*alpha**4-8*alpha**2+1

    rho = sp.symbols("rho", real=True)
    block = sp.zeros(4)
    for row, column, value in (
        (0, 1, 1), (0, 2, 1), (0, 3, -rho),
        (1, 3, 1), (2, 3, 1),
    ):
        block[row, column] = value
        block[column, row] = -value
    q_direction = sp.diag(block, -block)
    k_matrix = sp.simplify(system["delta"] / sp.I)
    tangent = specialize(
        sp.expand(system["x3"]*q_direction*system["x3"]
                  + k_matrix*q_direction*k_matrix/4)
    )
    tangent_equations: list[sp.Expr] = []
    for row in range(8):
        for column in range(row+1, 8):
            expression = sp.expand(tangent[row, column])
            if expression == 0:
                continue
            if any(
                sp.expand(expression-other) == 0
                or sp.expand(expression+other) == 0
                for other in tangent_equations
            ):
                continue
            tangent_equations.append(expression)
    if len(tangent_equations) != 6:
        raise AssertionError(
            f"P485 requires six tangent orbits, found {len(tangent_equations)}"
        )

    return {
        "system": system,
        "alpha": alpha,
        "rho": rho,
        "sine_substitution": sine_substitution,
        "N0": n0,
        "X_plus": x_plus,
        "X_minus": x_minus,
        "C": cross,
        "active_equations": active_equations,
        "branch_equation": branch_equation,
        "alpha_polynomial": alpha_polynomial,
        "tangent_equations": tangent_equations,
        "Q_direction": q_direction,
    }


def certificate_box(context: dict[str, Any]) -> dict[sp.Symbol, Interval]:
    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
    radius = Fraction(certificate["radius"])
    center = [Fraction(value) for value in certificate["center"]]
    result = {
        symbol: (value-radius, value+radius)
        for symbol, value in zip(context["system"]["variables"], center)
    }
    alpha_interval = tuple(
        Fraction(value) for value in certificate["sine_boxes"][0]
    )
    result[context["alpha"]] = alpha_interval  # type: ignore[assignment]
    return result


def center_substitution(context: dict[str, Any]) -> dict[sp.Symbol, float]:
    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
    center = [float(Fraction(value)) for value in certificate["center"]]
    result = dict(zip(context["system"]["variables"], center))
    result[context["alpha"]] = float(np.sin(np.pi/8))
    return result


def iv_add(left: Interval, right: Interval) -> Interval:
    return left[0]+right[0], left[1]+right[1]


def iv_mul(left: Interval, right: Interval) -> Interval:
    products = (
        left[0]*right[0], left[0]*right[1],
        left[1]*right[0], left[1]*right[1],
    )
    return min(products), max(products)


def iv_pow(value: Interval, exponent: int) -> Interval:
    if exponent == 0:
        return Fraction(1), Fraction(1)
    result: Interval = (Fraction(1), Fraction(1))
    base = value
    power = exponent
    while power:
        if power & 1:
            result = iv_mul(result, base)
        base = iv_mul(base, base)
        power //= 2
    return result


def polynomial_interval(
    expression: sp.Expr,
    symbols: tuple[sp.Symbol, ...],
    box: dict[sp.Symbol, Interval],
) -> Interval:
    polynomial = sp.Poly(sp.expand(expression), *symbols, domain=sp.QQ)
    total: Interval = (Fraction(0), Fraction(0))
    for exponents, coefficient in polynomial.terms():
        coefficient_fraction = Fraction(int(coefficient.p), int(coefficient.q))
        term: Interval = (coefficient_fraction, coefficient_fraction)
        for symbol, exponent in zip(symbols, exponents):
            if exponent:
                term = iv_mul(term, iv_pow(box[symbol], int(exponent)))
        total = iv_add(total, term)
    return total


def best_active_pivot(context: dict[str, Any]) -> dict[str, Any]:
    """Return the best floating-conditioned exact 3x3 A/B/u pivot."""

    A, B, u = context["system"]["variables"][:3]
    equations = context["active_equations"]
    coefficient_matrix = sp.Matrix([
        [sp.diff(value, variable) for variable in (A, B, u)]
        for value in equations
    ])
    substitution = center_substitution(context)
    numeric = np.asarray(coefficient_matrix.subs(substitution), dtype=float)
    candidates: list[tuple[float, tuple[int, int, int]]] = []
    for rows in itertools.combinations(range(12), 3):
        candidates.append((
            abs(float(np.linalg.det(numeric[list(rows), :]))), rows
        ))
    score, rows = max(candidates)
    pivot_matrix = coefficient_matrix[list(rows), :]
    determinant = sp.expand(pivot_matrix.det())
    return {
        "rows": rows,
        "floating_absolute_determinant": score,
        "coefficient_matrix": coefficient_matrix,
        "pivot_matrix": pivot_matrix,
        "determinant": determinant,
    }


def reduced_branch_system(context: dict[str, Any]) -> dict[str, Any]:
    """Eliminate A/B/u exactly and select the positive standard branch."""

    system = context["system"]
    A, B, u, L, a, b, c, d, e, f, g, h, i = system["variables"]
    alpha = context["alpha"]
    equations = context["active_equations"]
    pivot = best_active_pivot(context)
    rows = pivot["rows"]
    coefficient_matrix = pivot["coefficient_matrix"]
    pivot_matrix = pivot["pivot_matrix"]
    determinant = pivot["determinant"]
    zero_normalizer = {A: 0, B: 0, u: 0}
    constants = sp.Matrix([sp.expand(value.subs(zero_normalizer)) for value in equations])
    pivot_constants = constants[list(rows), :]
    adjugate = pivot_matrix.adjugate()
    normalizer_numerators = sp.Matrix([
        sp.expand(value) for value in (-adjugate*pivot_constants)
    ])
    remaining_rows = [index for index in range(12) if index not in rows]
    residual_numerators = []
    for index in remaining_rows:
        row = coefficient_matrix[index, :]
        numerator = sp.expand(
            determinant*constants[index] + (row*normalizer_numerators)[0]
        )
        residual_numerators.append(numerator)

    branch_substitution = {b: L-sp.Rational(36, 125)*alpha}
    reduced_equations = [
        sp.expand(value.subs(branch_substitution))
        for value in residual_numerators
    ]
    reduced_determinant = sp.expand(determinant.subs(branch_substitution))
    reduced_normalizer_numerators = [
        sp.expand(value.subs(branch_substitution))
        for value in normalizer_numerators
    ]
    reduced_variables = (a, c, d, e, f, g, h, i, alpha, L)
    primitive_equations = []
    for expression in reduced_equations:
        polynomial = sp.Poly(expression, *reduced_variables, domain=sp.QQ)
        _, primitive = polynomial.primitive()
        primitive_equations.append(primitive.as_expr())
    return {
        "pivot": pivot,
        "remaining_rows": remaining_rows,
        "normalizer_numerators": reduced_normalizer_numerators,
        "pivot_determinant": reduced_determinant,
        "branch_substitution": branch_substitution,
        "variables": reduced_variables,
        "equations": primitive_equations,
        "alpha_polynomial": context["alpha_polynomial"],
    }


def tangent_consistency(context: dict[str, Any]) -> dict[str, Any]:
    """Build exact consistency minors for the six affine rho equations."""

    rho = context["rho"]
    tangent = context["tangent_equations"]
    coefficients = [sp.expand(sp.diff(value, rho)) for value in tangent]
    constants = [sp.expand(value.subs(rho, 0)) for value in tangent]
    substitution = center_substitution(context)
    magnitudes = [abs(float(value.subs(substitution))) for value in coefficients]
    reference = max(range(len(coefficients)), key=magnitudes.__getitem__)
    consistency = []
    for index in range(len(tangent)):
        if index == reference:
            continue
        consistency.append(sp.expand(
            coefficients[reference]*constants[index]
            - coefficients[index]*constants[reference]
        ))
    return {
        "reference": reference,
        "coefficients": coefficients,
        "constants": constants,
        "consistency_polynomials": consistency,
        "floating_coefficient_magnitudes": magnitudes,
    }


def polynomial_stats(
    expressions: Iterable[sp.Expr], symbols: tuple[sp.Symbol, ...]
) -> list[dict[str, int]]:
    result = []
    for expression in expressions:
        polynomial = sp.Poly(sp.expand(expression), *symbols, domain=sp.QQ)
        result.append({
            "total_degree": int(polynomial.total_degree()),
            "term_count": len(polynomial.terms()),
        })
    return result
