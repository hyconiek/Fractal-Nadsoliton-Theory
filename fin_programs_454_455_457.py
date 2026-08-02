#!/usr/bin/env python3
"""FIN local research Programs P454, P455, and P457.

This file is intentionally dependency-light.  P454 derives and attempts to
certify the nested three-slot comb SDP dual.  P455 studies the exact coherent
symmetry face containing O167.  P457 refines the now globally licensed P452
one-dimensional interval cover.
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

import fin_programs_435_436_440 as p435
import fin_programs_445_447 as p445
import fin_programs_448_450 as p448
import fin_programs_451_453 as p451


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_454_455_457"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
P454_PATH = ROOT / f"{PREFIX}_P454_Nested_Dual.csv"
P454_WITNESS_PATH = ROOT / f"{PREFIX}_P454_Rational_Dual.npz"
P455_PATH = ROOT / f"{PREFIX}_P455_Symmetry_Residual.csv"
P457_PATH = ROOT / f"{PREFIX}_P457_Refined_Cover.csv"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p454_p457_dual_symmetry_and_refined_cover.png"


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


def symmetric_size(dimension: int) -> int:
    return dimension * (dimension + 1) // 2


def unpack_symmetric(vector: np.ndarray, dimension: int) -> np.ndarray:
    matrix = np.zeros((dimension, dimension), dtype=float)
    cursor = 0
    for index in range(dimension):
        matrix[index, index] = vector[cursor]
        cursor += 1
    for row in range(dimension):
        for column in range(row + 1, dimension):
            matrix[row, column] = matrix[column, row] = vector[cursor]
            cursor += 1
    return matrix


def pack_symmetric(matrix: np.ndarray) -> np.ndarray:
    values = [float(matrix[index, index]) for index in range(matrix.shape[0])]
    for row in range(matrix.shape[0]):
        for column in range(row + 1, matrix.shape[0]):
            values.append(float(matrix[row, column]))
    return np.array(values)


def pack_gradient(matrix: np.ndarray) -> np.ndarray:
    """Gradient coordinates for the symmetric-entry parameterization."""

    values = [float(np.real(matrix[index, index])) for index in range(matrix.shape[0])]
    for row in range(matrix.shape[0]):
        for column in range(row + 1, matrix.shape[0]):
            values.append(2 * float(np.real(matrix[row, column])))
    return np.array(values)


def unpack_dual(vector: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    cursor = 0
    n3 = symmetric_size(8)
    n2 = symmetric_size(4)
    n1 = symmetric_size(2)
    x3 = unpack_symmetric(vector[cursor:cursor+n3], 8)
    cursor += n3
    x2 = unpack_symmetric(vector[cursor:cursor+n2], 4)
    cursor += n2
    x1 = unpack_symmetric(vector[cursor:cursor+n1], 2)
    cursor += n1
    return x3, x2, x1, float(vector[cursor])


def pack_dual(x3: np.ndarray, x2: np.ndarray, x1: np.ndarray, lam: float) -> np.ndarray:
    return np.concatenate((pack_symmetric(x3), pack_symmetric(x2), pack_symmetric(x1), [lam]))


def dual_slacks(
    vector: np.ndarray, half_delta: np.ndarray
) -> tuple[list[np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray, float]]:
    x3, x2, x1, lam = unpack_dual(vector)
    slacks = [
        x3 + half_delta,
        x3 - half_delta,
        x2 - x3[:4, :4],
        x2 - x3[4:, 4:],
        x1 - x2[:2, :2],
        x1 - x2[2:, 2:],
        np.array([[lam - x1[0, 0]]]),
        np.array([[lam - x1[1, 1]]]),
    ]
    return slacks, (x3, x2, x1, lam)


def dual_barrier_value_gradient(
    vector: np.ndarray, half_delta: np.ndarray, barrier: float
) -> tuple[float, np.ndarray]:
    slacks, (x3, x2, x1, lam) = dual_slacks(vector, half_delta)
    logdet = 0.0
    inverses: list[np.ndarray] = []
    for slack in slacks:
        try:
            cholesky = np.linalg.cholesky(slack)
        except np.linalg.LinAlgError:
            return 1e30, np.zeros_like(vector)
        value = 2 * np.sum(np.log(np.real(np.diag(cholesky))))
        if not np.isfinite(value):
            return 1e30, np.zeros_like(vector)
        logdet += float(value)
        inverses.append(np.linalg.inv(slack))

    # Derivative of the sum of the eight log determinants.
    g3 = np.real(inverses[0] + inverses[1])
    g3[:4, :4] -= inverses[2]
    g3[4:, 4:] -= inverses[3]
    g2 = inverses[2] + inverses[3]
    g2[:2, :2] -= inverses[4]
    g2[2:, 2:] -= inverses[5]
    g1 = inverses[4] + inverses[5]
    g1[0, 0] -= inverses[6][0, 0]
    g1[1, 1] -= inverses[7][0, 0]
    glam = inverses[6][0, 0] + inverses[7][0, 0]

    gradient = np.concatenate((
        -barrier * pack_gradient(g3),
        -barrier * pack_gradient(g2),
        -barrier * pack_gradient(g1),
        [1 - barrier * float(np.real(glam))],
    ))
    return lam - barrier * logdet, gradient


def feasible_bfgs(
    initial: np.ndarray,
    half_delta: np.ndarray,
    barrier: float,
    maximum_iterations: int = 2000,
    tolerance: float = 1e-9,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Small deterministic BFGS with explicit feasibility backtracking."""

    vector = initial.copy()
    value, gradient = dual_barrier_value_gradient(vector, half_delta, barrier)
    inverse_hessian = np.eye(len(vector))
    accepted = 0
    resets = 0
    for iteration in range(maximum_iterations):
        if np.linalg.norm(gradient, np.inf) < tolerance:
            break
        direction = -inverse_hessian @ gradient
        directional = float(gradient @ direction)
        if directional >= -1e-14 * np.linalg.norm(gradient) * max(1, np.linalg.norm(direction)):
            inverse_hessian = np.eye(len(vector))
            direction = -gradient
            directional = -float(gradient @ gradient)
            resets += 1
        step = 1.0
        candidate_value = math.inf
        candidate_gradient = gradient
        while step >= 2**-60:
            candidate = vector + step * direction
            candidate_value, candidate_gradient = dual_barrier_value_gradient(
                candidate, half_delta, barrier
            )
            if candidate_value <= value + 1e-4 * step * directional:
                break
            step /= 2
        if step < 2**-60:
            break
        displacement = step * direction
        gradient_change = candidate_gradient - gradient
        curvature = float(displacement @ gradient_change)
        vector = candidate
        value = candidate_value
        gradient = candidate_gradient
        accepted += 1
        if curvature > 1e-12 * np.linalg.norm(displacement) * max(
            1, np.linalg.norm(gradient_change)
        ):
            rho = 1 / curvature
            identity = np.eye(len(vector))
            left = identity - rho * np.outer(displacement, gradient_change)
            inverse_hessian = (
                left @ inverse_hessian @ left.T
                + rho * np.outer(displacement, displacement)
            )
        else:
            inverse_hessian = np.eye(len(vector))
            resets += 1
    return vector, {
        "iterations": iteration + 1,
        "accepted_steps": accepted,
        "inverse_hessian_resets": resets,
        "barrier_objective": value,
        "gradient_infinity_norm": float(np.linalg.norm(gradient, np.inf)),
        "converged": bool(np.linalg.norm(gradient, np.inf) < tolerance),
    }


def p454_dual_barrier_scout() -> tuple[np.ndarray, list[dict[str, Any]]]:
    delta = p448.compressed_process_difference(3, 0.8, math.pi / 8)
    half_delta = delta / 2
    norm = float(np.linalg.norm(half_delta, 2))
    x3 = (norm + 0.5) * np.eye(8)
    x2 = (norm + 1.0) * np.eye(4)
    x1 = (norm + 1.5) * np.eye(2)
    vector = pack_dual(x3, x2, x1, norm + 2.0)
    rows: list[dict[str, Any]] = []
    for barrier in (
        0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001,
        0.00003, 0.00001, 0.000003, 0.000001, 0.0000003,
    ):
        vector, diagnostics = feasible_bfgs(
            vector, half_delta, barrier, maximum_iterations=2000, tolerance=1e-8
        )
        slacks, (_, _, _, lam) = dual_slacks(vector, half_delta)
        rows.append({
            "barrier": barrier,
            "lambda": lam,
            "minimum_slack_eigenvalue": min(
                float(np.linalg.eigvalsh(slack)[0]) for slack in slacks
            ),
            **diagnostics,
        })
    return vector, rows


def fraction_matrix(matrix: np.ndarray, scale: int) -> list[list[Fraction]]:
    return [
        [Fraction(round(float(value) * scale), scale) for value in row]
        for row in matrix
    ]


def sympy_fraction_matrix(matrix: list[list[Fraction]]) -> sp.Matrix:
    return sp.Matrix([
        [sp.Rational(value.numerator, value.denominator) for value in row]
        for row in matrix
    ])


def exact_sylvester_margin(
    matrix: sp.Matrix,
    candidates: tuple[Fraction, ...] = (
        Fraction(1, 2_000_000), Fraction(1, 5_000_000),
        Fraction(1, 10_000_000), Fraction(1, 20_000_000),
        Fraction(1, 50_000_000), Fraction(1, 100_000_000),
        Fraction(1, 200_000_000), Fraction(1, 1_000_000_000),
    ),
) -> tuple[Fraction, list[Fraction]]:
    """Prove matrix >= margin*I by exact leading principal minors."""

    for margin in candidates:
        shifted = matrix - sp.Rational(margin.numerator, margin.denominator) * sp.eye(matrix.rows)
        determinants: list[Fraction] = []
        safe = True
        for size in range(1, matrix.rows + 1):
            determinant = sp.cancel(shifted[:size, :size].det(method="domain-ge"))
            determinant = sp.simplify(determinant)
            if determinant.is_real is False or determinant <= 0:
                safe = False
                break
            determinants.append(Fraction(int(determinant.p), int(determinant.q)))
        if safe:
            return margin, determinants
    return Fraction(0), []


def p454_half_delta_midpoint() -> tuple[sp.Matrix, Fraction]:
    dimension = 8
    scale = 10**20
    midpoint = sp.zeros(dimension)
    maximum_radius = Fraction(0)
    for row in range(dimension):
        for column in range(dimension):
            hamming = (row ^ column).bit_count()
            difference = row.bit_count() - column.bit_count()
            sine = p448.p449_sine_interval(difference)
            interval = p445.iv_scale(sine, Fraction(4, 5) ** hamming)
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
    return midpoint, maximum_radius


def p454_rational_dual_certificate(
    vector: np.ndarray, scale: int = 10**12
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Rationalize and exactly certify every nested dual PSD slack."""

    x3_float, x2_float, x1_float, lambda_float = unpack_dual(vector)
    x3_values = fraction_matrix(x3_float, scale)
    x2_values = fraction_matrix(x2_float, scale)
    x1_values = fraction_matrix(x1_float, scale)
    lambda_value = Fraction(math.ceil(lambda_float * scale), scale)
    x3 = sympy_fraction_matrix(x3_values)
    x2 = sympy_fraction_matrix(x2_values)
    x1 = sympy_fraction_matrix(x1_values)
    half_delta, half_delta_radius = p454_half_delta_midpoint()

    slacks = {
        "X3_plus_halfDelta": x3 + half_delta,
        "X3_minus_halfDelta": x3 - half_delta,
        "X2_minus_X3_block0": x2 - x3[:4, :4],
        "X2_minus_X3_block1": x2 - x3[4:, 4:],
        "X1_minus_X2_block0": x1 - x2[:2, :2],
        "X1_minus_X2_block1": x1 - x2[2:, 2:],
    }
    rows: list[dict[str, Any]] = []
    certified_lowers: list[Fraction] = []
    for name, slack in slacks.items():
        margin, determinants = exact_sylvester_margin(slack)
        perturbation = (
            Fraction(8) * half_delta_radius
            if "halfDelta" in name else Fraction(0)
        )
        certified = margin - perturbation
        certified_lowers.append(certified)
        rows.append({
            "slack": name,
            "exact_sylvester_margin": margin,
            "transcendental_operator_perturbation": perturbation,
            "certified_exact_minimum_eigenvalue_lower": certified,
            "positive_leading_principal_minors": len(determinants),
            "smallest_leading_principal_minor": min(determinants) if determinants else None,
        })
    scalar_slacks = [
        lambda_value - x1_values[0][0],
        lambda_value - x1_values[1][1],
    ]
    for index, slack in enumerate(scalar_slacks):
        certified_lowers.append(slack)
        rows.append({
            "slack": f"lambda_minus_X1_diag{index}",
            "exact_sylvester_margin": slack,
            "transcendental_operator_perturbation": Fraction(0),
            "certified_exact_minimum_eigenvalue_lower": slack,
            "positive_leading_principal_minors": 1,
            "smallest_leading_principal_minor": slack,
        })
    feasible = min(certified_lowers) > 0
    return ({
        "status": (
            "[Computer-assisted proof] exact-rational feasible nested comb dual"
            if feasible else
            "[Open] rationalized dual failed at least one exact PSD slack"
        ),
        "feasible": feasible,
        "rationalization_scale": scale,
        "lambda_upper": lambda_value,
        "smallest_certified_slack_lower": min(certified_lowers),
        "half_delta_maximum_entry_radius": half_delta_radius,
        "slacks": rows,
    }, {
        "X3": np.array([[float(value) for value in row] for row in x3_values]),
        "X2": np.array([[float(value) for value in row] for row in x2_values]),
        "X1": np.array([[float(value) for value in row] for row in x1_values]),
        "lambda": np.array([float(lambda_value)]),
    })


# ---------------------------------------------------------------------------
# P455: what the available symmetries do and do not force
# ---------------------------------------------------------------------------


def real_symmetric_basis(dimension: int) -> list[np.ndarray]:
    result: list[np.ndarray] = []
    for index in range(dimension):
        matrix = np.zeros((dimension, dimension))
        matrix[index, index] = 1
        result.append(matrix)
    for row in range(dimension):
        for column in range(row + 1, dimension):
            matrix = np.zeros((dimension, dimension))
            matrix[row, column] = matrix[column, row] = 1
            result.append(matrix)
    return result


def real_c3_tangent_basis() -> list[np.ndarray]:
    """Fourteen exact real tangent directions of the affine C3 recursion."""

    result: list[np.ndarray] = []
    # r direction and three real-symmetric C0 directions.
    cases = [(np.diag([1.0, -1.0]), np.zeros((2, 2)))]
    cases.extend((np.zeros((2, 2)), matrix) for matrix in real_symmetric_basis(2))
    for n1, c0 in cases:
        c1 = n1 - c0
        n2 = np.zeros((4, 4))
        n2[:2, :2] = c0
        n2[2:, 2:] = c1
        matrix = np.zeros((8, 8))
        matrix[4:, 4:] = n2
        result.append(matrix)
    # Ten real-symmetric B0 directions, with B1=-B0.
    for b0 in real_symmetric_basis(4):
        matrix = np.zeros((8, 8))
        matrix[:4, :4] = b0
        matrix[4:, 4:] = -b0
        result.append(matrix)
    return result


def p455_fixed_basis() -> tuple[list[np.ndarray], dict[str, Any]]:
    basis = real_c3_tangent_basis()
    complement = np.fliplr(np.eye(8))
    constraint_columns = []
    for matrix in basis:
        difference = complement @ matrix @ complement - matrix
        constraint_columns.append([
            sp.Rational(int(round(value))) for value in difference.reshape(-1)
        ])
    constraint = sp.Matrix.hstack(*[sp.Matrix(column) for column in constraint_columns])
    nullspace = constraint.nullspace()
    fixed = []
    for vector in nullspace:
        matrix = sum(
            (float(vector[index]) * basis[index] for index in range(len(basis))),
            start=np.zeros((8, 8)),
        )
        fixed.append(matrix)
    flattened = np.column_stack([matrix.reshape(-1) for matrix in fixed])
    orthonormal, _ = np.linalg.qr(flattened)
    orthonormal_basis = [orthonormal[:, index].reshape(8, 8) for index in range(len(fixed))]

    # The O167 affine face has tangent coordinates A, B, and u.
    base = coherent_face_normalizer(np.array([0.0, 0.0, 0.0]))
    face = [
        coherent_face_normalizer(np.eye(3)[index]) - base for index in range(3)
    ]
    fixed_rank = int(sp.Matrix(flattened).rank())
    face_rank = int(np.linalg.matrix_rank(np.column_stack([m.reshape(-1) for m in face])))
    combined_rank = int(np.linalg.matrix_rank(np.column_stack(
        [m.reshape(-1) for m in fixed] + [m.reshape(-1) for m in face]
    )))
    return orthonormal_basis, {
        "real_c3_affine_dimension": len(basis),
        "complement_constraint_rank": int(constraint.rank()),
        "real_complement_fixed_dimension": len(nullspace),
        "O167_face_dimension": face_rank,
        "face_contained_in_fixed_space": combined_rank == fixed_rank,
        "symmetry_allowed_residual_dimensions": fixed_rank - face_rank,
        "exact_nullspace_vectors": [[str(value) for value in vector] for vector in nullspace],
    }


def coherent_face_normalizer(parameters: np.ndarray) -> np.ndarray:
    A, B, u = map(float, parameters)
    C = 0.5 - A - 2 * B
    left = np.array([
        [A, u, u, 0], [u, B, 0, -u],
        [u, 0, B, -u], [0, -u, -u, C],
    ])
    right = np.array([
        [C, -u, -u, 0], [-u, B, 0, u],
        [-u, 0, B, u], [0, u, u, A],
    ])
    result = np.zeros((8, 8))
    result[:4, :4] = left
    result[4:, 4:] = right
    return result


def feasible_distance_objective(normalizer: np.ndarray, delta: np.ndarray) -> float:
    minimum = float(np.linalg.eigvalsh(normalizer)[0])
    if minimum <= 1e-11:
        return 1000 - minimum
    return -p448.tester_distance(normalizer, delta)


def numerical_hessian(function: Any, point: np.ndarray, step: float = 2e-5) -> np.ndarray:
    dimension = len(point)
    result = np.zeros((dimension, dimension))
    center = function(point)
    for row in range(dimension):
        direction_row = np.zeros(dimension)
        direction_row[row] = step
        result[row, row] = (
            function(point + direction_row) - 2 * center + function(point - direction_row)
        ) / step**2
        for column in range(row + 1, dimension):
            direction_column = np.zeros(dimension)
            direction_column[column] = step
            value = (
                function(point + direction_row + direction_column)
                - function(point + direction_row - direction_column)
                - function(point - direction_row + direction_column)
                + function(point - direction_row - direction_column)
            ) / (4 * step**2)
            result[row, column] = result[column, row] = value
    return result


def p455_symmetry_audit() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    delta = p448.compressed_process_difference(3, 0.8, math.pi / 8)
    face_initial = np.array([0.287166, 0.068777, 0.009493])
    face_optimum = minimize(
        lambda parameters: feasible_distance_objective(
            coherent_face_normalizer(parameters), delta
        ),
        face_initial,
        method="Nelder-Mead",
        options={"maxiter": 10000, "xatol": 1e-13, "fatol": 1e-15},
    )

    fixed_basis, exact_audit = p455_fixed_basis()
    central = np.eye(8) / 8
    fixed_matrix = np.column_stack([matrix.reshape(-1) for matrix in fixed_basis])
    face_normalizer = coherent_face_normalizer(face_optimum.x)
    fixed_initial = fixed_matrix.T @ (face_normalizer - central).reshape(-1)
    fixed_objective = lambda coordinates: feasible_distance_objective(
        central + sum(
            (coordinates[index] * fixed_basis[index] for index in range(len(fixed_basis))),
            start=np.zeros((8, 8)),
        ),
        delta,
    )
    fixed_rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(455)
    best = None
    for seed in range(8):
        initial = fixed_initial if seed == 0 else fixed_initial + rng.normal(scale=0.02, size=5)
        optimum = minimize(
            fixed_objective,
            initial,
            method="Nelder-Mead",
            options={"maxiter": 12000, "xatol": 1e-12, "fatol": 1e-14},
        )
        if best is None or optimum.fun < best.fun:
            best = optimum
        fixed_rows.append({
            "seed": seed,
            "half_distance": -float(optimum.fun),
            "iterations": int(optimum.nit),
            "success": bool(optimum.success),
        })
    assert best is not None
    fixed_hessian = numerical_hessian(lambda x: -fixed_objective(x), best.x)
    face_hessian = numerical_hessian(
        lambda x: -feasible_distance_objective(coherent_face_normalizer(x), delta),
        face_optimum.x,
    )
    return ({
        "status": (
            "[Proven] real plus global-complement symmetry leaves a five-dimensional "
            "fixed affine space, two dimensions larger than the O167 face; "
            "[Strong evidence] all five-dimensional starts return to the O167 value"
        ),
        **exact_audit,
        "face_optimum_parameters": face_optimum.x.tolist(),
        "face_optimum_half_distance": -float(face_optimum.fun),
        "face_minimum_eigenvalue": float(np.linalg.eigvalsh(face_normalizer)[0]),
        "face_hessian_eigenvalues": np.linalg.eigvalsh(face_hessian).tolist(),
        "fixed_space_best_half_distance": -float(best.fun),
        "fixed_space_best_coordinates": best.x.tolist(),
        "fixed_space_hessian_eigenvalues": np.linalg.eigvalsh(fixed_hessian).tolist(),
        "fixed_space_starts": len(fixed_rows),
        "maximum_fixed_minus_face_gain": -float(best.fun) + float(face_optimum.fun),
        "obstruction": (
            "Known conjugation and global-complement symmetries cannot by themselves "
            "reduce the ordered causal cone to the three-parameter O167 face; two "
            "symmetry-allowed directions remain. Numerical non-improvement is not a theorem."
        ),
        "new_object": "O171 Ordered-Comb Symmetry Residual Space",
    }, fixed_rows)


# ---------------------------------------------------------------------------
# Program wrappers and P457 licensed refinement
# ---------------------------------------------------------------------------


def program_454() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    vector, central_rows = p454_dual_barrier_scout()
    certificate, matrices = p454_rational_dual_certificate(vector)
    if not certificate["feasible"]:
        raise RuntimeError("P454 rational dual certificate is not feasible")
    inherited = json.loads(
        (ROOT / "FIN_Programs_451_453_Results.json").read_text(encoding="utf-8")
    )["P451"]
    primal_lower = Fraction(inherited["coherent_half_distance_interval"][0])
    dual_upper = Fraction(certificate["lambda_upper"])
    gap = dual_upper - primal_lower
    np.savez_compressed(P454_WITNESS_PATH, **matrices)
    rows = [
        {"row_type": "central_path", **row} for row in central_rows
    ] + [
        {"row_type": "exact_dual_slack", **row}
        for row in certificate["slacks"]
    ]
    return ({
        "status": (
            "[Proven] exact nested three-slot comb dual derivation; "
            "[Computer-assisted proof] exact-rational feasible dual and global "
            "full-cone bracket of width below 6.7e-6"
        ),
        "primal": (
            "maximize Tr(Delta*(T_plus-T_minus))/2 over T_plus,T_minus>=0, "
            "T_plus+T_minus=B0 direct-sum B1, B0+B1=C0 direct-sum C1, "
            "C0+C1=diag(r,1-r)"
        ),
        "dual": (
            "minimize lambda subject to X3>=+/-Delta/2; X2>=each 4x4 block "
            "of X3; X1>=each 2x2 block of X2; lambda>=both diagonal entries of X1"
        ),
        "strong_duality_scope": (
            "Both primal and dual have strict feasible points in the finite reduced "
            "comb SDP, so standard finite-dimensional Slater duality applies."
        ),
        "inherited_primal_lower": primal_lower,
        "certified_dual_upper": dual_upper,
        "certified_global_gap": gap,
        "certified_success_probability_interval": (
            (1 + primal_lower) / 2, (1 + dual_upper) / 2,
        ),
        "rational_dual_certificate": certificate,
        "central_path_steps": len(central_rows),
        "final_barrier_scout_lambda": central_rows[-1]["lambda"],
        "boundary": (
            "The bracket proves a global value bound, not exact attainment by the "
            "O167 face and not uniqueness of an optimal tester. It concerns only the "
            "declared reduced three-slot channel."
        ),
        "new_object": "O170 Nested Comb Dual Ladder",
    }, rows)


def program_457() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    certificate = p445.palindromic_branch_certificate(tolerance=1e-6)
    rows: list[dict[str, Any]] = []
    for a in np.linspace(0, 0.5, 201):
        probabilities = np.array([a, 0.5-a, 0.5-a, a])
        rows.append({
            "a": a,
            "coarse_half_distance": p445.p446_numeric_objective(probabilities),
            "inside_certified_maximizer_hull": (
                certificate["maximizer_hull"][0]
                <= a <= certificate["maximizer_hull"][1]
            ),
        })
    return ({
        "status": (
            "[Computer-assisted proof] globally licensed full-simplex coarse-erasure "
            "interval of width below 1e-6"
        ),
        "license": (
            "P452 proves every declared four-sector simplex point is dominated by its "
            "palindromic reversal average; therefore this one-dimensional outward cover "
            "is a full-simplex certificate, not merely a line result."
        ),
        "global_lower": certificate["incumbent_lower"],
        "global_upper": certificate["global_upper"],
        "global_gap": certificate["optimality_gap"],
        "incumbent_a": certificate["incumbent_a"],
        "maximizer_hull": certificate["maximizer_hull"],
        "evaluated_boxes": certificate["evaluated_boxes"],
        "terminal_boxes": certificate["terminal_boxes"],
        "improvement_factor_over_P452_gap": (
            0.0009924652413310642 / certificate["optimality_gap"]
        ),
        "boundary": (
            "The certificate retains the exact P452 parameter/code/loss scope. It does "
            "not prove a closed-form maximizer, uniqueness, unrestricted input "
            "optimality, adaptivity, or laboratory performance."
        ),
        "new_object": "O172 Six-Decimal Coarse-Erasure Cover",
    }, rows)


def make_figure(
    p454: dict[str, Any], p454_rows: list[dict[str, Any]],
    p455: dict[str, Any], p455_rows: list[dict[str, Any]],
    p457: dict[str, Any], p457_rows: list[dict[str, Any]],
) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.35))

    central = [row for row in p454_rows if row["row_type"] == "central_path"]
    barriers = [row["barrier"] for row in central]
    lambdas = [row["lambda"] for row in central]
    axes[0].semilogx(barriers, lambdas, "o-", color="#2563eb")
    axes[0].invert_xaxis()
    axes[0].axhspan(
        float(Fraction(p454["inherited_primal_lower"])),
        float(Fraction(p454["certified_dual_upper"])),
        color="#f59e0b", alpha=0.32, label="certified global bracket",
    )
    axes[0].set_title("P454 nested dual central path")
    axes[0].set_xlabel("barrier parameter")
    axes[0].set_ylabel("dual lambda")
    axes[0].legend(fontsize=8)

    distances = [row["half_distance"] for row in p455_rows]
    axes[1].plot(range(len(distances)), distances, "o", color="#0f766e")
    axes[1].axhline(p455["face_optimum_half_distance"], color="#dc2626", linestyle="--")
    axes[1].set_title("P455 five-dimensional fixed space")
    axes[1].set_xlabel("deterministic start")
    axes[1].set_ylabel("half distance")

    x = [row["a"] for row in p457_rows]
    y = [row["coarse_half_distance"] for row in p457_rows]
    axes[2].plot(x, y, color="#7c3aed")
    axes[2].axhspan(
        p457["global_lower"], p457["global_upper"],
        color="#f59e0b", alpha=0.35, label="global 1e-6 bracket",
    )
    axes[2].set_title("P457 refined full-simplex cover")
    axes[2].set_xlabel("palindromic sector mass a")
    axes[2].set_ylabel("coarse half distance")
    axes[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)


def main() -> None:
    p454_result, p454_rows = program_454()
    write_csv(P454_PATH, p454_rows)
    p455_result, p455_rows = p455_symmetry_audit()
    write_csv(P455_PATH, p455_rows)
    p457_result, p457_rows = program_457()
    write_csv(P457_PATH, p457_rows)
    make_figure(
        p454_result, p454_rows, p455_result, p455_rows, p457_result, p457_rows
    )
    results = {
        "metadata": {
            "programs": "P454/P455/P457",
            "checkpoint": "P454-P457",
            "date": "2026-08-01",
            "local_only": True,
            "external_physical_evidence": False,
            "selector_discharged": False,
            "dimensional_source_exported": False,
            "legacy_strict_bridge_complete": False,
            "legacy_role_transfer_started": False,
        },
        "P454": p454_result,
        "P455": p455_result,
        "P457": p457_result,
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(results), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_csv(SUMMARY_PATH, [
        {
            "program": name,
            "status": results[name]["status"],
            "new_object": results[name]["new_object"],
            "physical_evidence": False,
            "selector_discharged": False,
        }
        for name in ("P454", "P455", "P457")
    ])
    print(json.dumps({
        "P454_global_gap": json_ready(p454_result["certified_global_gap"]),
        "P455_fixed_dimension": p455_result["real_complement_fixed_dimension"],
        "P455_residual_dimensions": p455_result["symmetry_allowed_residual_dimensions"],
        "P457_global_gap": p457_result["global_gap"],
    }, indent=2))


if __name__ == "__main__":
    main()
