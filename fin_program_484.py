#!/usr/bin/env python3
"""FIN local Program P484: adversarial audit of the complex phase-face claim.

The proof removes the numerical-nullspace caveat of P474.  It uses the exact
swap/complement representation of the P471 matrices, the positive Riccati
identity proved in P473, and the odd-dimensional fixed-axis theorem for an
orientation-preserving metric-orthogonal 3x3 map.

The finite computations below prove the representation formulas and audit a
high-precision candidate tangent.  They also expose the remaining exact
causality/intersection obligation, so floating residuals are not promoted to
an exact phase-face theorem.
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
import sympy as sp

import fin_programs_448_450 as p448
import fin_programs_471_472_473 as p471


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Program_484"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
ROWS_PATH = ROOT / f"{PREFIX}_Symmetry_Phase_Face.csv"
WITNESS_PATH = ROOT / f"{PREFIX}_Phase_Face_Witness.npz"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p484_exact_complex_phase_face.png"


def json_ready(value: Any) -> Any:
    if isinstance(value, sp.MatrixBase):
        return json_ready(value.tolist())
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, sp.Basic):
        return str(value)
    if isinstance(value, np.generic):
        return json_ready(value.item())
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
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
            writer.writerow({key: json_ready(row.get(key, "")) for key in keys})


def six_dimensional_symmetry_basis() -> sp.Matrix:
    """Swap-symmetric, complement-even/odd basis, ordered +++---."""

    root_two = sp.sqrt(2)
    basis = sp.zeros(8, 6)
    for row, value in ((0, 1/root_two), (7, 1/root_two)):
        basis[row, 0] = value
    for row in (1, 2, 5, 6):
        basis[row, 1] = sp.Rational(1, 2)
    for row, value in ((3, 1/root_two), (4, 1/root_two)):
        basis[row, 2] = value
    for row, value in ((0, 1/root_two), (7, -1/root_two)):
        basis[row, 3] = value
    for row, value in (
        (1, sp.Rational(1, 2)), (2, sp.Rational(1, 2)),
        (5, -sp.Rational(1, 2)), (6, -sp.Rational(1, 2)),
    ):
        basis[row, 4] = value
    for row, value in ((3, 1/root_two), (4, -1/root_two)):
        basis[row, 5] = value
    return basis


def exact_representation_audit(system: dict[str, Any]) -> dict[str, Any]:
    """Prove the exact block forms used by the P484 analytic argument."""

    transform = six_dimensional_symmetry_basis()
    normalizer = system["normalizer"]
    x3 = system["x3"]
    delta = system["delta"]
    k_matrix = sp.simplify(delta / sp.I)
    reduced_n = sp.simplify(transform.T * normalizer * transform)
    reduced_x = sp.simplify(transform.T * x3 * transform)
    reduced_k = sp.simplify(transform.T * k_matrix * transform)
    n_plus = reduced_n[:3, :3]
    n_minus = reduced_n[3:, 3:]
    x_plus = reduced_x[:3, :3]
    x_minus = reduced_x[3:, 3:]
    c_block = reduced_k[:3, 3:]

    # The swap-antisymmetric complement-parity vectors are not columns of the
    # six-dimensional basis.  They carry the remaining two one-dimensional
    # blocks and expose the scalar standard-sector equation.
    eye = sp.eye(8)
    antisymmetric_plus = (
        eye[:, 1] - eye[:, 2] + eye[:, 6] - eye[:, 5]
    ) / 2
    antisymmetric_minus = (
        eye[:, 1] - eye[:, 2] - eye[:, 6] + eye[:, 5]
    ) / 2
    standard_x_plus = sp.simplify((antisymmetric_plus.T*x3*antisymmetric_plus)[0])
    standard_x_minus = sp.simplify((antisymmetric_minus.T*x3*antisymmetric_minus)[0])
    standard_k = sp.simplify((antisymmetric_plus.T*k_matrix*antisymmetric_minus)[0])

    exact_checks = {
        "basis_orthonormal": sp.simplify(transform.T*transform-sp.eye(6)) == sp.zeros(6),
        "normalizer_parity_offblock_zero": (
            reduced_n[:3, 3:] == sp.zeros(3)
            and reduced_n[3:, :3] == sp.zeros(3)
        ),
        "normalizer_plus_equals_minus": sp.simplify(n_plus-n_minus) == sp.zeros(3),
        "X_parity_offblock_zero": (
            reduced_x[:3, 3:] == sp.zeros(3)
            and reduced_x[3:, :3] == sp.zeros(3)
        ),
        "K_parity_diagonal_blocks_zero": (
            reduced_k[:3, :3] == sp.zeros(3)
            and reduced_k[3:, 3:] == sp.zeros(3)
        ),
        "K_lower_cross_is_negative_transpose": (
            sp.simplify(reduced_k[3:, :3] + c_block.T) == sp.zeros(3)
        ),
        "standard_X_blocks_equal": sp.simplify(standard_x_plus-standard_x_minus) == 0,
        "standard_X_value": sp.simplify(standard_x_plus-(system["variables"][3]-system["variables"][5])) == 0,
        "standard_K_value": sp.simplify(standard_k + sp.Rational(72,125)*system["sines"][0]) == 0,
    }
    if not all(exact_checks.values()):
        raise AssertionError(f"failed exact P484 representation check: {exact_checks}")
    return {
        "checks": exact_checks,
        "N_plus": n_plus,
        "X_plus": x_plus,
        "X_minus": x_minus,
        "C": c_block,
        "standard_X": standard_x_plus,
        "standard_K": standard_k,
    }


def odd_dimension_causality_counterexample() -> dict[str, Any]:
    """Exact counterexample to an over-broad odd-dimensional inference.

    Active Riccati congruence and an orientation-preserving metric rotation do
    not alone force the fixed axis into the equal-endpoint causal plane.
    Additional shared-entry/standard-sector identities are indispensable.
    """

    identity = sp.eye(3)
    axis = sp.Matrix([1, 0, 2])
    rotation = 2 * (axis * axis.T) / (axis.dot(axis)) - identity
    cross = 2 * rotation.T
    checks = {
        "rotation_orthogonal": sp.simplify(rotation*rotation.T-identity) == sp.zeros(3),
        "rotation_orientation_preserving": sp.det(rotation) == 1,
        "fixed_axis": sp.simplify(rotation.T*axis-axis) == sp.zeros(3, 1),
        "axis_not_causal_equal_endpoint": axis[0] != axis[2],
        "plus_Riccati_block": sp.simplify(identity-cross*cross.T/4) == sp.zeros(3),
        "minus_Riccati_block": sp.simplify(identity-cross.T*cross/4) == sp.zeros(3),
    }
    if not all(checks.values()):
        raise AssertionError(f"failed P484 exact counterexample: {checks}")
    return {
        "checks": checks,
        "N0": identity,
        "X_plus": identity,
        "X_minus": identity,
        "C": cross,
        "B": rotation,
        "fixed_axis": axis,
    }
def high_precision_root_and_axis(
    system: dict[str, Any], precision: int = 100
) -> tuple[list[sp.Expr], list[sp.Expr], sp.Expr]:
    """Locate the exact branch at high precision; never use it as a zero proof."""

    variables = system["variables"]
    sines = system["sines"]
    s1 = sp.sqrt(2-sp.sqrt(2))/2
    exact_sines = (s1, sp.sqrt(2)/2, sp.sqrt(2+sp.sqrt(2))/2)
    equations = [value.subs(dict(zip(sines, exact_sines))) for value in system["equations"]]
    initial = list(np.load(
        ROOT / "FIN_Programs_471_472_473_P473_Root_Box.npz"
    )["center"])
    root = list(sp.nsolve(
        equations,
        variables,
        initial,
        tol=sp.Float("1e-90"),
        maxsteps=200,
        prec=precision,
        solver="mdnewton",
    ))

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
    tangent = sp.expand(
        system["x3"]*q_direction*system["x3"] + k_matrix*q_direction*k_matrix/4
    )
    unique: list[sp.Expr] = []
    for row in range(8):
        for column in range(row+1, 8):
            expression = sp.expand(tangent[row, column])
            if expression == 0:
                continue
            if any(
                sp.expand(expression-other) == 0 or sp.expand(expression+other) == 0
                for other in unique
            ):
                continue
            unique.append(expression)
    substitutions = dict(zip(variables, root))
    substitutions.update(dict(zip(sines, exact_sines)))
    ratios: list[sp.Expr] = []
    for expression in unique:
        polynomial = sp.Poly(expression, rho)
        coefficient = polynomial.coeff_monomial(rho).subs(substitutions)
        constant = polynomial.coeff_monomial(1).subs(substitutions)
        ratios.append(sp.N(-constant/coefficient, precision-20))
    spread = sp.N(max(ratios)-min(ratios), precision-20)
    return root, ratios, spread


def numeric_phase_face(
    system: dict[str, Any], point: np.ndarray
) -> dict[str, Any]:
    """Audit the 3x3 metric rotation and its causal fixed axis."""

    normalizer, x3 = p471.structured_numeric_matrices(point)
    delta = p448.compressed_process_difference(3, 0.8, math.pi/8)
    k_matrix = np.real(delta/1j)
    transform = np.asarray(six_dimensional_symmetry_basis(), dtype=float)
    reduced_n = transform.T @ normalizer @ transform
    reduced_x = transform.T @ x3 @ transform
    reduced_k = transform.T @ k_matrix @ transform
    n0 = reduced_n[:3, :3]
    x_plus = reduced_x[:3, :3]
    x_minus = reduced_x[3:, 3:]
    c_block = reduced_k[:3, 3:]
    a_map = np.linalg.solve(x_plus, c_block)/2
    b_map = np.linalg.solve(x_minus, c_block.T)/2
    values, vectors = np.linalg.eig(b_map.T)
    axis = np.real(vectors[:, np.argmin(np.abs(values-1))])
    axis /= axis[0]
    skew_plus = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0],
    ])
    parity_q = np.zeros((6, 6))
    parity_q[:3, :3] = skew_plus
    parity_q[3:, 3:] = skew_plus
    q_direction = transform @ parity_q @ transform.T
    tangent = x3@q_direction@x3 + k_matrix@q_direction@k_matrix/4
    inherited_lower = json.loads(
        (ROOT / "FIN_Programs_471_472_473_Results.json").read_text(encoding="utf-8")
    )["P473"]["normalizer_box_positive_lower"]
    positive_lower = Fraction(inherited_lower)
    operator_norm = float(np.linalg.norm(q_direction, 2))
    safe_t = float(positive_lower)/(2*operator_norm)
    return {
        "A_times_B_residual": float(np.linalg.norm(a_map@b_map-np.eye(3))),
        "metric_orthogonality_residual": float(np.linalg.norm(b_map@n0@b_map.T-n0)),
        "B_determinant": float(np.linalg.det(b_map)),
        "B_eigenvalues": np.linalg.eigvals(b_map),
        "fixed_axis": axis,
        "axis_balance_residual": float(abs(axis[0]-axis[2])),
        "tangent_residual_infinity_norm": float(np.linalg.norm(tangent, np.inf)),
        "Q_operator_norm": operator_norm,
        "inherited_exact_normalizer_lower": positive_lower,
        "displayed_safe_half_interval": safe_t,
        "Q": q_direction,
        "N": normalizer,
        "X3": x3,
        "Delta": delta,
    }


def make_figure(audit: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    normalizer = audit["N"]
    q_direction = audit["Q"]
    delta = audit["Delta"]
    x3 = audit["X3"]
    safe = audit["displayed_safe_half_interval"]
    parameters = np.linspace(-safe, safe, 61)
    eigenvalues = []
    riccati = []
    values = []
    for parameter in parameters:
        candidate = normalizer + 1j*parameter*q_direction
        eigenvalues.append(float(np.min(np.linalg.eigvalsh(candidate))))
        riccati.append(float(np.linalg.norm(
            x3@candidate@x3-delta@candidate@delta/4, np.inf
        )))
        values.append(p448.tester_distance(candidate, delta))
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].plot(parameters, eigenvalues, color="tab:blue")
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set_xlabel("phase-face parameter t")
    axes[0].set_ylabel("minimum eigenvalue of N(t)")
    axes[0].set_title("P484 numerically admissible phase segment")
    axes[0].grid(alpha=0.3)
    center = values[len(values)//2]
    axes[1].plot(parameters, np.asarray(values)-center, label="objective change")
    axes[1].plot(parameters, riccati, label="Riccati residual")
    axes[1].set_yscale("symlog", linthresh=1e-16)
    axes[1].set_xlabel("phase-face parameter t")
    axes[1].set_title("Numerical audit of conditional invariance")
    axes[1].grid(alpha=0.3)
    axes[1].legend()
    figure.tight_layout()
    figure.savefig(FIGURE_PATH, dpi=180)
    plt.close(figure)


def main() -> None:
    system = p471.polynomial_system()
    representation = exact_representation_audit(system)
    counterexample = odd_dimension_causality_counterexample()
    root, ratios, spread = high_precision_root_and_axis(system)
    point = np.asarray([float(value) for value in root])
    audit = numeric_phase_face(system, point)
    make_figure(audit)
    np.savez_compressed(
        WITNESS_PATH,
        point=point,
        Q=audit["Q"],
        N=audit["N"],
        X3=audit["X3"],
        Delta=audit["Delta"],
        fixed_axis=audit["fixed_axis"],
    )
    rows = [
        {"row_type": "exact_representation", "check": key, "passed": value}
        for key, value in representation["checks"].items()
    ]
    rows.extend({
        "row_type": "high_precision_rho",
        "equation_orbit": index,
        "rho": sp.N(value, 80),
    } for index, value in enumerate(ratios))
    rows.extend((
        {"row_type": "metric_rotation", "quantity": "AB_residual", "value": audit["A_times_B_residual"]},
        {"row_type": "metric_rotation", "quantity": "metric_orthogonality_residual", "value": audit["metric_orthogonality_residual"]},
        {"row_type": "metric_rotation", "quantity": "det_B", "value": audit["B_determinant"]},
        {"row_type": "causal_axis", "quantity": "axis_balance_residual", "value": audit["axis_balance_residual"]},
        {"row_type": "causal_axis", "quantity": "tangent_residual", "value": audit["tangent_residual_infinity_norm"]},
    ))
    write_csv(ROWS_PATH, rows)
    result = {
        "metadata": {
            "program": "P484",
            "execution_mode": "local analytical and computational research only",
            "network_used": False,
            "laboratory_data_used": False,
            "external_audit_used": False,
            "selector_boundary": "QW-2191 remains open.",
            "kernel_boundary": "No legacy/strict substitution or role transfer is made.",
            "physical_boundary": "No dimensional or experimentally validated physics is claimed.",
        },
        "status": (
            "[Strong evidence; exact proof obligation open] the P473 O167 "
            "optimizer has a nontrivial complex affine phase-face candidate; "
            "full-cone optimizer uniqueness is not yet rigorously refuted"
        ),
        "exact_representation_checks": representation["checks"],
        "exact_odd_dimension_counterexample": {
            "checks": counterexample["checks"],
            "B": counterexample["B"],
            "fixed_axis": counterexample["fixed_axis"],
            "meaning": (
                "With N0=X+=X-=I3 and C=2B^T, both active Riccati blocks hold "
                "and B is orientation-preserving orthogonal, yet its fixed axis "
                "(1,0,2) is outside the causal equal-endpoint plane. Thus odd "
                "dimension and active Riccati congruence alone are insufficient."
            ),
        },
        "conditional_analytic_theorem": (
            "Swap and complement reduce the active real sector to equal positive "
            "3x3 normalizer blocks N0, dual blocks X+/X-, and a cross block C. "
            "Positive Riccati congruence gives X-=C^T X+^-1 C/4. Hence "
            "B=X-^-1 C^T/2 is orientation-preserving N0-orthogonal. In odd "
            "dimension B has a fixed skew form. IF its fixed covector lies in the "
            "causal equal-endpoint plane, transforming back gives nonzero real skew Q with "
            "X3 Q X3+K Q K/4=0. Therefore N(t)=N+i t Q is causal, positive for "
            "all sufficiently small nonzero t, obeys the same Riccati equation, "
            "and attains the same exact dual value."
        ),
        "open_exact_step": (
            "The exact shared-entry and standard-sector equations have not been "
            "shown by a symbolic identity or interval existence certificate to "
            "force the fixed covector into the causal equal-endpoint plane. The "
            "six tangent equations agree to high precision, but this is not a zero proof."
        ),
        "full_cone_optimizer_uniqueness_refuted": False,
        "exact_nonzero_optimal_segment_proved": False,
        "conditional_nonzero_optimal_segment_theorem": True,
        "numerical_phase_face_dimension_lower_bound": 1,
        "complete_optimal_face_dimension_proved": False,
        "rho_high_precision_values": ratios,
        "rho_high_precision_spread": spread,
        "rho_spread_boundary": (
            "Zero here means equality after the displayed 80-digit rounding. "
            "A 180-digit rerun resolves differences of order 1e-102 caused by "
            "the finite nonlinear-solve tolerance; neither observation proves exact equality."
        ),
        "numeric_metric_rotation": {
            key: value for key, value in audit.items()
            if key not in {"Q", "N", "X3", "Delta"}
        },
        "conditional_positivity_statement": (
            "If the candidate Q is proved to be the exact causal tangent, then "
            "after scaling Q to operator norm one, the exact P473 lower bound "
            f"mu={audit['inherited_exact_normalizer_lower']} guarantees "
            "N+i t Q>0 for every |t|<mu."
        ),
        "interpretation": (
            "O167 selects a unique local real polynomial representative and an "
            "exact dual support; whether the complete complex operational normalizer "
            "is unique remains open, with strong evidence for nonuniqueness."
        ),
        "boundary": (
            "The audit does not prove nonuniqueness or an exact phase face. It "
            "isolates one exact missing identity/certificate. Even conditional on "
            "that identity, it would not classify the complete face, turn the phase "
            "into a selector or observable, supply physical units, or add evidence."
        ),
        "new_object": "O185 Conditional Odd-Dimensional Complex Comb Phase Face",
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(result), indent=2, ensure_ascii=False)+"\n",
        encoding="utf-8",
    )
    print(json.dumps(json_ready(result), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
