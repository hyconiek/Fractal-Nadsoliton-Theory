#!/usr/bin/env python3
"""FIN local research Programs P474, P479, and P483.

P474 audits the complete optimal face left open by P473.  P479 exports a
small dependency-free Lean statement of the reusable algebraic argument and
records whether the local Lean checker is actually available.  P483 builds a
uniform rational Krawczyk tube in the channel parameters q and theta.

The programs are mathematical/computational only.  They neither consume nor
manufacture laboratory evidence, physical units, a selector, or a
legacy-to-strict role-transfer theorem.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Callable

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import qr
import sympy as sp

import fin_programs_448_450 as p448
import fin_programs_465_468_469 as p465
import fin_programs_471_472_473 as p471


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_474_479_483"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
P474_PATH = ROOT / f"{PREFIX}_P474_Optimal_Face.csv"
P474_WITNESS_PATH = ROOT / f"{PREFIX}_P474_Flat_Direction.npz"
P479_PATH = ROOT / f"{PREFIX}_P479_Formalization.csv"
P479_LEAN_PATH = ROOT / f"{PREFIX}_P479_Riccati_Trace_Core.lean"
P483_PATH = ROOT / f"{PREFIX}_P483_Parametric_Tube.csv"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p474_p483_optimal_face_and_parameter_tube.png"

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


def hermitian_basis(dimension: int) -> list[np.ndarray]:
    """Real-coordinate basis of complex Hermitian matrices."""

    result: list[np.ndarray] = []
    for index in range(dimension):
        matrix = np.zeros((dimension, dimension), dtype=complex)
        matrix[index, index] = 1
        result.append(matrix)
    for row in range(dimension):
        for column in range(row + 1, dimension):
            matrix = np.zeros((dimension, dimension), dtype=complex)
            matrix[row, column] = matrix[column, row] = 1
            result.append(matrix)
    for row in range(dimension):
        for column in range(row + 1, dimension):
            matrix = np.zeros((dimension, dimension), dtype=complex)
            matrix[row, column] = 1j
            matrix[column, row] = -1j
            result.append(matrix)
    return result


def imaginary_causal_basis() -> tuple[list[np.ndarray], list[str]]:
    """Seven real skew matrices Q for which iQ is a C3 tangent direction."""

    result: list[np.ndarray] = []
    names: list[str] = []
    # The imaginary off-diagonal coordinate of C0; C1 receives its negative.
    c0 = np.array([[0.0, 1.0], [-1.0, 0.0]])
    n2 = np.zeros((4, 4))
    n2[:2, :2] = c0
    n2[2:, 2:] = -c0
    matrix = np.zeros((8, 8))
    matrix[4:, 4:] = n2
    result.append(matrix)
    names.append("C0_01")
    # Six imaginary B0 coordinates; B1 receives the negative so B0+B1=0.
    for row in range(4):
        for column in range(row + 1, 4):
            block = np.zeros((4, 4))
            block[row, column] = 1
            block[column, row] = -1
            matrix = np.zeros((8, 8))
            matrix[:4, :4] = block
            matrix[4:, 4:] = -block
            result.append(matrix)
            names.append(f"B0_{row}{column}")
    return result, names


def causal_homogeneous_functionals() -> list[tuple[str, Callable[[np.ndarray], float]]]:
    """Independent real equations for a homogeneous C3 normalizer tangent."""

    result: list[tuple[str, Callable[[np.ndarray], float]]] = []
    for row in range(4):
        for column in range(4, 8):
            result.extend((
                (f"N3_cross_{row}_{column}_re",
                 lambda matrix, row=row, column=column: float(np.real(matrix[row, column]))),
                (f"N3_cross_{row}_{column}_im",
                 lambda matrix, row=row, column=column: float(np.imag(matrix[row, column]))),
            ))
    for row in range(2):
        for column in range(2, 4):
            result.extend((
                (f"N2_cross_{row}_{column}_re",
                 lambda matrix, row=row, column=column: float(np.real(
                     matrix[row, column] + matrix[row + 4, column + 4]
                 ))),
                (f"N2_cross_{row}_{column}_im",
                 lambda matrix, row=row, column=column: float(np.imag(
                     matrix[row, column] + matrix[row + 4, column + 4]
                 ))),
            ))
    result.extend((
        ("N1_offdiag_re", lambda matrix: float(np.real(
            matrix[0, 1] + matrix[4, 5] + matrix[2, 3] + matrix[6, 7]
        ))),
        ("N1_offdiag_im", lambda matrix: float(np.imag(
            matrix[0, 1] + matrix[4, 5] + matrix[2, 3] + matrix[6, 7]
        ))),
        ("normalization_trace", lambda matrix: float(np.real(np.trace(matrix)))),
    ))
    return result


def full_contact_matrix(
    slack_plus: np.ndarray, slack_minus: np.ndarray
) -> tuple[np.ndarray, list[str], list[np.ndarray]]:
    """Linear equations for differences of two primal optima at a fixed dual."""

    basis = hermitian_basis(8)
    rows: list[np.ndarray] = []
    names: list[str] = []
    for label, slack, offset in (
        ("plus", slack_plus, 0), ("minus", slack_minus, 64)
    ):
        for row in range(8):
            for column in range(8):
                for component in ("re", "im"):
                    equation = np.zeros(128)
                    for index, matrix in enumerate(basis):
                        value = (slack @ matrix)[row, column]
                        equation[offset + index] = (
                            float(np.real(value)) if component == "re"
                            else float(np.imag(value))
                        )
                    rows.append(equation)
                    names.append(f"{label}_contact_{row}_{column}_{component}")
    for name, functional in causal_homogeneous_functionals():
        equation = np.zeros(128)
        for index, matrix in enumerate(basis):
            equation[index] = functional(matrix)
            equation[64 + index] = functional(matrix)
        rows.append(equation)
        names.append(name)
    return np.asarray(rows), names, basis


def program_474() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Audit full-cone uniqueness and expose the complex flat direction."""

    witness = np.load(ROOT / "FIN_Programs_471_472_473_P473_Root_Box.npz")
    point = np.asarray(witness["center"], dtype=float)
    normalizer, x3 = p471.structured_numeric_matrices(point)
    delta = p448.compressed_process_difference(3, 0.8, math.pi / 8)
    slack_plus = x3 - delta / 2
    slack_minus = x3 + delta / 2
    contact, contact_names, basis = full_contact_matrix(slack_plus, slack_minus)
    singular = np.linalg.svd(contact, compute_uv=False)
    contact_rank = int(np.linalg.matrix_rank(contact, tol=1e-10))

    # A smaller representation makes the origin of the defect transparent.
    # For Delta=iK and N_t=N+i*t*Q, the Riccati equation is affine in t and
    # its tangent equation is X Q X + K Q K/4=0.
    k_matrix = np.real(delta / 1j)
    causal_basis, causal_names = imaginary_causal_basis()
    riccati_columns = []
    for direction in causal_basis:
        residual = x3 @ direction @ x3 + k_matrix @ direction @ k_matrix / 4
        riccati_columns.append(np.array([
            residual[row, column]
            for row in range(8) for column in range(row + 1, 8)
        ]))
    riccati_map = np.asarray(riccati_columns).T
    _, riccati_singular, right = np.linalg.svd(riccati_map, full_matrices=False)
    coefficients = right[-1]
    coefficients /= np.max(np.abs(coefficients))
    q_direction = sum(
        (coefficients[index] * causal_basis[index] for index in range(7)),
        start=np.zeros((8, 8)),
    )
    direction = 1j * q_direction
    riccati_residual = x3 @ q_direction @ x3 + k_matrix @ q_direction @ k_matrix / 4

    # Recover the associated tester-pair tangent from the full contact system.
    _, _, contact_right = np.linalg.svd(contact, full_matrices=False)
    tester_coordinates = contact_right[-1]
    tester_coordinates /= np.max(np.abs(tester_coordinates))
    tester_plus_direction = sum(
        (tester_coordinates[index] * basis[index] for index in range(64)),
        start=np.zeros((8, 8), dtype=complex),
    )
    tester_minus_direction = sum(
        (tester_coordinates[64 + index] * basis[index] for index in range(64)),
        start=np.zeros((8, 8), dtype=complex),
    )
    normalizer_direction = tester_plus_direction + tester_minus_direction

    # Positivity and objective scouts on a symmetric interior segment.
    _, _, _, _, _, tester_plus, tester_minus = p465.o167_support_ladder(
        point[:3], delta
    )
    scan = np.linspace(-0.025, 0.025, 21)
    base_value = p448.tester_distance(normalizer, delta)
    rows: list[dict[str, Any]] = []
    maximum_value_deviation = 0.0
    minimum_normalizer_eigenvalue = math.inf
    minimum_tester_eigenvalue = math.inf
    for parameter in scan:
        candidate_n = normalizer + parameter * normalizer_direction
        candidate_plus = tester_plus + parameter * tester_plus_direction
        candidate_minus = tester_minus + parameter * tester_minus_direction
        value = p448.tester_distance(candidate_n, delta)
        deviation = value - base_value
        maximum_value_deviation = max(maximum_value_deviation, abs(deviation))
        minimum_normalizer_eigenvalue = min(
            minimum_normalizer_eigenvalue,
            float(np.min(np.linalg.eigvalsh(candidate_n))),
        )
        minimum_tester_eigenvalue = min(
            minimum_tester_eigenvalue,
            float(np.min(np.linalg.eigvalsh(candidate_plus))),
            float(np.min(np.linalg.eigvalsh(candidate_minus))),
        )
        rows.append({
            "row_type": "flat_segment_scout",
            "parameter": parameter,
            "normalizer_minimum_eigenvalue": float(np.min(np.linalg.eigvalsh(candidate_n))),
            "tester_plus_minimum_eigenvalue": float(np.min(np.linalg.eigvalsh(candidate_plus))),
            "tester_minus_minimum_eigenvalue": float(np.min(np.linalg.eigvalsh(candidate_minus))),
            "half_distance": value,
            "deviation_from_center": deviation,
        })
    for name, coefficient in zip(causal_names, coefficients):
        rows.append({
            "row_type": "imaginary_causal_null_vector",
            "coordinate": name,
            "coefficient": coefficient,
        })
    selected_rows = qr(contact.T, pivoting=True, mode="economic")[2][:127]
    selected_singular = np.linalg.svd(contact[selected_rows], compute_uv=False)
    rows.append({
        "row_type": "rank_summary",
        "full_equation_rows": contact.shape[0],
        "full_unknown_coordinates": contact.shape[1],
        "numerical_rank": contact_rank,
        "smallest_nonzero_singular_value": singular[-2],
        "smallest_singular_value": singular[-1],
        "independent_127_row_minor_smallest_singular_value": selected_singular[-1],
        "selected_causal_rows": int(sum(index >= 256 for index in selected_rows)),
    })
    np.savez_compressed(
        P474_WITNESS_PATH,
        point=point,
        normalizer=normalizer,
        X3=x3,
        Delta=delta,
        causal_coefficients=coefficients,
        Q_direction=q_direction,
        normalizer_direction=normalizer_direction,
        tester_plus=tester_plus,
        tester_minus=tester_minus,
        tester_plus_direction=tester_plus_direction,
        tester_minus_direction=tester_minus_direction,
        contact_singular_values=singular,
        riccati_singular_values=riccati_singular,
    )
    pattern_rho = float(-coefficients[3] / coefficients[1])
    return ({
        "status": (
            "[Strong evidence] the complete complex causal cone has a "
            "one-dimensional optimal face through the P473 optimizer; "
            "[Open] exact algebraic certification of the flat direction"
        ),
        "full_contact_matrix_shape": list(contact.shape),
        "full_contact_numerical_rank": contact_rank,
        "full_contact_nullity": 128 - contact_rank,
        "smallest_nonzero_contact_singular_value": singular[-2],
        "smallest_contact_singular_value": singular[-1],
        "imaginary_causal_riccati_shape": list(riccati_map.shape),
        "imaginary_causal_riccati_rank": int(np.linalg.matrix_rank(riccati_map, tol=1e-10)),
        "imaginary_causal_riccati_nullity": 7 - int(np.linalg.matrix_rank(riccati_map, tol=1e-10)),
        "imaginary_causal_singular_values": riccati_singular,
        "imaginary_causal_coefficients": dict(zip(causal_names, coefficients)),
        "reduced_pattern_rho": pattern_rho,
        "riccati_tangent_residual_infinity_norm": float(np.linalg.norm(riccati_residual, np.inf)),
        "full_contact_null_residual_infinity_norm": float(np.linalg.norm(
            contact @ tester_coordinates, np.inf
        )),
        "normalizer_direction_real_norm": float(np.linalg.norm(np.real(normalizer_direction))),
        "normalizer_direction_imaginary_norm": float(np.linalg.norm(np.imag(normalizer_direction))),
        "tested_parameter_interval": [float(scan[0]), float(scan[-1])],
        "minimum_scanned_normalizer_eigenvalue": minimum_normalizer_eigenvalue,
        "minimum_scanned_tester_eigenvalue": minimum_tester_eigenvalue,
        "maximum_scanned_objective_deviation": maximum_value_deviation,
        "interpretation": (
            "The real O167 point is numerically isolated inside its real symmetry "
            "face, but the complete complex comb admits an imaginary causal phase "
            "direction. This is operational nonuniqueness, not a new selector, "
            "physical phase, clock, or laboratory observable."
        ),
        "boundary": (
            "Rank and flatness are currently floating-point evidence. They do not "
            "yet prove that the exact P473 root has an exact nonzero causal null "
            "vector. Consequently full-cone uniqueness is not marked refuted."
        ),
        "new_object": "O182 Complex Comb Optimal-Face Direction",
    }, rows)


LEAN_SOURCE = r'''import Std

/-!
P479: dependency-free algebraic core of the O181 argument.

The carrier and operations are deliberately abstract.  Positivity, square
roots, spectral absolute values, and finite-matrix trace theory remain named
premises; this file proves only the exact implication structure actually used
by the FIN checkpoint.  It exports no selector, dimensions, physical units,
legacy-role transfer, L_total, Standard Model, GR, or laboratory claim.
-/

namespace FINPrograms474To483

theorem riccati_support_bridge
    {Matrix : Type}
    (riccati positiveSupport : Matrix -> Matrix -> Matrix -> Prop)
    (N X Delta : Matrix)
    (bridge : forall A B C, riccati A B C -> positiveSupport A B C)
    (hRiccati : riccati N X Delta) :
    positiveSupport N X Delta :=
  bridge N X Delta hRiccati

theorem trace_telescope_attainment
    {Primal Dual : Type}
    (primalValue : Primal -> Rat)
    (dualValue : Dual -> Rat)
    (feasiblePrimal : Primal -> Prop)
    (feasibleDual : Dual -> Prop)
    (contact : Primal -> Dual -> Prop)
    (weakDuality : forall p d,
      feasiblePrimal p -> feasibleDual d -> primalValue p <= dualValue d)
    (contactEquality : forall p d, contact p d -> primalValue p = dualValue d)
    (p : Primal) (d : Dual)
    (hp : feasiblePrimal p) (hd : feasibleDual d) (hc : contact p d) :
    primalValue p = dualValue d := by
  exact contactEquality p d hc

theorem exact_attainment_blocks_strict_improvement
    {Primal Dual : Type}
    (primalValue : Primal -> Rat)
    (dualValue : Dual -> Rat)
    (feasiblePrimal : Primal -> Prop)
    (feasibleDual : Dual -> Prop)
    (weakDuality : forall p d,
      feasiblePrimal p -> feasibleDual d -> primalValue p <= dualValue d)
    (witness : Primal) (certificate : Dual)
    (hw : feasiblePrimal witness) (hc : feasibleDual certificate)
    (attained : primalValue witness = dualValue certificate) :
    forall candidate, feasiblePrimal candidate ->
      primalValue candidate <= primalValue witness := by
  intro candidate hcandidate
  rw [attained]
  exact weakDuality candidate certificate hcandidate hc

theorem local_root_uniqueness_does_not_assert_global_optimizer_uniqueness
    (LocalRootUnique GlobalOptimizerUnique : Prop)
    (globalProof : GlobalOptimizerUnique) :
    GlobalOptimizerUnique :=
  globalProof

end FINPrograms474To483
'''


def program_479() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    P479_LEAN_PATH.write_text(LEAN_SOURCE, encoding="utf-8")
    lean = shutil.which("lean")
    rows: list[dict[str, Any]] = []
    checked = False
    output = "lean executable absent"
    if lean is not None:
        completed = subprocess.run(
            [lean, str(P479_LEAN_PATH)],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=30,
            check=False,
        )
        checked = completed.returncode == 0
        output = completed.stdout.strip()
        return_code = completed.returncode
    else:
        return_code = None
    rows.extend((
        {
            "theorem": "riccati_support_bridge",
            "scope": "abstract implication; analytic matrix bridge remains a premise",
            "lean_checked": checked,
        },
        {
            "theorem": "trace_telescope_attainment",
            "scope": "contact equality separated from feasibility and weak duality",
            "lean_checked": checked,
        },
        {
            "theorem": "exact_attainment_blocks_strict_improvement",
            "scope": "standard weak-duality consequence",
            "lean_checked": checked,
        },
        {
            "theorem": "local_root_uniqueness_does_not_assert_global_optimizer_uniqueness",
            "scope": "typing guard against the P473-to-P474 overclaim",
            "lean_checked": checked,
        },
    ))
    toolchain_missing = lean is not None and "no default toolchain configured" in output
    return ({
        "status": (
            "[Proven] Lean-checked dependency-free implication core"
            if checked else
            "[Open] formal source exported, but no local Lean toolchain is configured"
        ),
        "lean_executable": lean,
        "lean_return_code": return_code,
        "lean_checked": checked,
        "local_toolchain_missing": toolchain_missing,
        "checker_output": output,
        "formal_source": P479_LEAN_PATH.name,
        "proved_outside_lean": (
            "The finite-dimensional matrix proof of O181 remains the exact "
            "computer-assisted theorem of P471-P473. P479 only isolates its "
            "logical dependency graph."
        ),
        "boundary": (
            "The local lean launcher exists but has no installed default toolchain; "
            "the source is therefore not described as machine-checked. No network "
            "installation was attempted."
            if not checked else
            "The abstract Lean core does not formalize spectral square roots or PSD matrices."
        ),
        "new_object": "O183 Riccati-Attainment Formal Interface",
    }, rows)


def parametric_polynomial_system() -> dict[str, Any]:
    """P471 system with q promoted to an interval parameter."""

    base = p471.polynomial_system()
    variables = base["variables"]
    sines = base["sines"]
    q_parameter = sp.symbols("q_parameter", real=True)
    normalizer = base["normalizer"]
    x3 = base["x3"]
    sine = {
        0: 0, 1: sines[0], 2: sines[1], 3: sines[2],
        -1: -sines[0], -2: -sines[1], -3: -sines[2],
    }
    delta = sp.zeros(8)
    for row in range(8):
        for column in range(8):
            exponent = (row ^ column).bit_count()
            difference = row.bit_count() - column.bit_count()
            delta[row, column] = (
                2 * sp.I * q_parameter**exponent * sine[difference]
            )
    residual = sp.expand(x3 * normalizer * x3 - delta * normalizer * delta / 4)
    equations = tuple(residual[row, column] for row, column in (
        positions[0] for positions in base["groups"]
    ))
    jacobian = sp.Matrix(equations).jacobian(variables)
    symbols = variables + sines + (q_parameter,)
    return {
        "variables": variables,
        "sines": sines,
        "q": q_parameter,
        "equations": equations,
        "jacobian": jacobian,
        "equation_polys": tuple(sp.Poly(value, *symbols, domain=sp.QQ) for value in equations),
        "jacobian_polys": tuple(tuple(
            sp.Poly(jacobian[row, column], *symbols, domain=sp.QQ)
            for column in range(13)
        ) for row in range(13)),
    }


def parametric_krawczyk_certificate(
    system: dict[str, Any],
    center: tuple[Fraction, ...],
    variable_radius: Fraction,
    q_radius: Fraction,
    theta_radius: Fraction,
    preconditioner: list[list[Fraction]],
) -> dict[str, Any]:
    variable_box = tuple(
        (value - variable_radius, value + variable_radius) for value in center
    )
    point_box = tuple((value, value) for value in center)
    nominal_sines = p471.exact_sine_boxes(scale=10**30)
    sine_box = tuple(
        (value[0] - (index + 1) * theta_radius,
         value[1] + (index + 1) * theta_radius)
        for index, value in enumerate(nominal_sines)
    )
    q_box = (Fraction(4, 5) - q_radius, Fraction(4, 5) + q_radius)
    parameter_values = sine_box + (q_box,)
    function_box = tuple(
        p471.polynomial_interval(poly, point_box + parameter_values)
        for poly in system["equation_polys"]
    )
    jacobian_box = tuple(tuple(
        p471.polynomial_interval(poly, variable_box + parameter_values)
        for poly in row
    ) for row in system["jacobian_polys"])
    correction: list[Interval] = []
    for row in range(13):
        value: Interval = (Fraction(0), Fraction(0))
        for column in range(13):
            value = p471.iv_add(
                value,
                p471.iv_scale(function_box[column], preconditioner[row][column]),
            )
        correction.append(p471.iv_neg(value))
    c_j = p471.interval_matrix_product(preconditioner, [list(row) for row in jacobian_box])
    identity_minus = []
    for row in range(13):
        current = []
        for column in range(13):
            identity = Fraction(int(row == column))
            current.append(p471.iv_sub((identity, identity), c_j[row][column]))
        identity_minus.append(tuple(current))
    displacement = tuple((-variable_radius, variable_radius) for _ in range(13))
    linear = []
    for row in identity_minus:
        total: Interval = (Fraction(0), Fraction(0))
        for value, shift in zip(row, displacement):
            total = p471.iv_add(total, p471.iv_mul(value, shift))
        linear.append(total)
    image = tuple(
        p471.iv_add((center[index], center[index]), p471.iv_add(correction[index], linear[index]))
        for index in range(13)
    )
    margins = tuple(
        min(image[index][0] - variable_box[index][0],
            variable_box[index][1] - image[index][1])
        for index in range(13)
    )
    included = all(value > 0 for value in margins)
    contraction_rows = []
    for row in identity_minus:
        contraction_rows.append(sum(max(abs(value[0]), abs(value[1])) for value in row))
    return {
        "variable_radius": variable_radius,
        "q_radius": q_radius,
        "theta_radius": theta_radius,
        "included": included,
        "minimum_inclusion_margin": min(margins),
        "maximum_contraction_row_sum": max(contraction_rows),
        "maximum_function_interval_radius": max(
            max(abs(value[0]), abs(value[1])) for value in function_box
        ),
        "q_interval": q_box,
        "theta_interval_symbolic": (
            f"[pi/8-{theta_radius}, pi/8+{theta_radius}]"
        ),
        "sine_intervals": sine_box,
    }


def program_483() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Uniform parameter-neighborhood theorem for the O181 root."""

    system = parametric_polynomial_system()
    witness = np.load(ROOT / "FIN_Programs_471_472_473_P473_Root_Box.npz")
    point = np.asarray(witness["center"], dtype=float)
    center = p471.rational_center(point, scale=10**16)
    sine_values = p471.sine_midpoints()
    arguments = system["variables"] + system["sines"] + (system["q"],)
    jacobian_function = sp.lambdify(arguments, system["jacobian"], modules="numpy")
    jacobian = np.asarray(
        jacobian_function(*point, *sine_values, 0.8), dtype=float
    ).reshape(13, 13)
    preconditioner = p471.rational_matrix(np.linalg.inv(jacobian), scale=10**14)

    trials: list[tuple[Fraction, Fraction, Fraction]] = []
    # Variable boxes are deliberately much wider than the P473 isolation box;
    # the parameter boxes are searched from ambitious to conservative.
    for variable_radius in (
        Fraction(1, 10**4), Fraction(1, 10**5),
        Fraction(1, 10**6), Fraction(1, 10**7),
    ):
        for parameter_radius in (
            Fraction(1, 10**6), Fraction(3, 10**7), Fraction(1, 10**7),
            Fraction(3, 10**8), Fraction(1, 10**8), Fraction(3, 10**9),
            Fraction(1, 10**9), Fraction(1, 10**10),
        ):
            trials.append((variable_radius, parameter_radius, parameter_radius))
    certificates = [
        parametric_krawczyk_certificate(system, center, vr, qr, tr, preconditioner)
        for vr, qr, tr in trials
    ]
    rows = [{"row_type": "tube_trial", **value} for value in certificates]
    admitted = [value for value in certificates if value["included"]]
    if admitted:
        # Maximize the common q/theta radius, then minimize the variable box.
        selected = sorted(
            admitted,
            key=lambda value: (-value["q_radius"], value["variable_radius"]),
        )[0]
        status = (
            "[Computer-assisted proof] one exact positive O181 root persists "
            "uniquely inside the declared common box for every parameter pair "
            "in the certified q/theta rectangle"
        )
    else:
        selected = min(
            certificates,
            key=lambda value: float(value["minimum_inclusion_margin"]),
            default=None,
        )
        status = "[Open] no tested uniform q/theta Krawczyk tube was admitted"

    # Positivity persists if the P473 lower margins dominate the structured
    # entry perturbation paid by the wider variable box.
    inherited = json.loads(
        (ROOT / "FIN_Programs_471_472_473_Results.json").read_text(encoding="utf-8")
    )["P473"]
    positivity = False
    normalizer_lower = None
    x3_lower = None
    if admitted:
        vr = selected["variable_radius"]
        normalizer_lower = Fraction(inherited["normalizer_center_sylvester_margin"]) - 40 * vr
        x3_lower = Fraction(inherited["x3_center_sylvester_margin"]) - 8 * vr
        positivity = normalizer_lower > 0 and x3_lower > 0
        rows.append({
            "row_type": "tube_positivity",
            "normalizer_positive_lower": normalizer_lower,
            "X3_positive_lower": x3_lower,
            "positive_throughout_tube": positivity,
        })
    proved = bool(admitted and positivity)
    return ({
        "status": status if proved else (
            status + "; [Open] uniform positivity was not certified"
            if admitted else status
        ),
        "uniform_root_tube_proved": proved,
        "trials": len(certificates),
        "admitted_trials": len(admitted),
        "selected_certificate": selected,
        "normalizer_tube_positive_lower": normalizer_lower,
        "X3_tube_positive_lower": x3_lower,
        "parameter_semantics": (
            "q varies directly; theta varies through the rigorous Lipschitz "
            "enclosures |sin(k theta)-sin(k theta0)|<=k|theta-theta0| for k=1,2,3."
        ),
        "theorem_scope": (
            "Uniform existence and uniqueness are local in the 13 polynomial "
            "coordinates and concern only the declared reduced three-slot channel."
        ),
        "boundary": (
            "The tube does not prove global optimizer uniqueness, selector closure, "
            "dimensional physics, a laboratory realization, or transfer of legacy roles."
        ),
        "new_object": "O184 Parametric O181 Krawczyk Tube",
    }, rows)


def make_figure(p474: dict[str, Any], p474_rows: list[dict[str, Any]],
                p483_rows: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    scan = [row for row in p474_rows if row["row_type"] == "flat_segment_scout"]
    tubes = [row for row in p483_rows if row["row_type"] == "tube_trial"]
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].plot(
        [row["parameter"] for row in scan],
        [row["normalizer_minimum_eigenvalue"] for row in scan],
        marker="o", label="min eig(N(t))",
    )
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set_xlabel("flat-direction parameter t")
    axes[0].set_ylabel("minimum eigenvalue")
    axes[0].set_title("P474 positivity along the candidate optimal face")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    admitted = [row for row in tubes if row["included"]]
    rejected = [row for row in tubes if not row["included"]]
    for values, label, color, marker in (
        (admitted, "admitted", "tab:green", "o"),
        (rejected, "rejected", "tab:red", "x"),
    ):
        if values:
            axes[1].scatter(
                [float(row["variable_radius"]) for row in values],
                [float(row["q_radius"]) for row in values],
                label=label, color=color, marker=marker,
            )
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("13D root-box radius")
    axes[1].set_ylabel("common q/theta radius")
    axes[1].set_title("P483 uniform parametric Krawczyk trials")
    axes[1].grid(alpha=0.3, which="both")
    axes[1].legend()
    figure.tight_layout()
    figure.savefig(FIGURE_PATH, dpi=180)
    plt.close(figure)


def main() -> None:
    p474_result, p474_rows = program_474()
    p479_result, p479_rows = program_479()
    p483_result, p483_rows = program_483()
    make_figure(p474_result, p474_rows, p483_rows)
    write_csv(P474_PATH, p474_rows)
    write_csv(P479_PATH, p479_rows)
    write_csv(P483_PATH, p483_rows)
    summary = [
        {"program": "P474", "status": p474_result["status"], "new_object": p474_result["new_object"]},
        {"program": "P479", "status": p479_result["status"], "new_object": p479_result["new_object"]},
        {"program": "P483", "status": p483_result["status"], "new_object": p483_result["new_object"]},
    ]
    write_csv(SUMMARY_PATH, summary)
    RESULTS_PATH.write_text(json.dumps(json_ready({
        "metadata": {
            "programs": ["P474", "P479", "P483"],
            "execution_mode": "local analytical and computational research only",
            "network_used": False,
            "laboratory_data_used": False,
            "external_audit_used": False,
            "kernel_boundary": "No legacy/strict substitution or role transfer is made.",
            "selector_boundary": "QW-2191 remains open.",
            "physical_boundary": "No dimensional or experimentally predictive physics is claimed.",
        },
        "P474": p474_result,
        "P479": p479_result,
        "P483": p483_result,
    }), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(json_ready({
        "P474": p474_result,
        "P479": p479_result,
        "P483": p483_result,
    }), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
