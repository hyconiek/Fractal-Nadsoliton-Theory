#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import numpy as np
from numpy.polynomial import Polynomial

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2569_s1519_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate import apd_rows

GEN = ROOT / "generated"
OUT = GEN / "p2576_s1526_apd_boundary_nullspace_discrete_sobolev_selector_dependence_audit.json"
MD = GEN / "p2576_s1526_apd_boundary_nullspace_discrete_sobolev_selector_dependence_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2574_APD_TWO_ENDPOINT_COMPATIBILITY": GEN / "p2574_s1524_apd_two_endpoint_boundary_compatibility_obstruction_audit.json",
    "P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE": GEN / "p2575_s1525_apd_augmented_boundary_nullspace_nonuniqueness_audit.json",
}
DOMAIN = list(range(12))
DERIVATIVE_ORDERS = [0, 1, 2, 3, 4, 12]
DISCRETE_GRID = [i / 2 for i in range(23)]
BOUNDARY_TARGETS = [
    {"name": "zero_zero_neumann", "left_slope_target": 0.0, "right_slope_target": 0.0},
    {"name": "equal_positive_slopes", "left_slope_target": 1.0e-6, "right_slope_target": 1.0e-6},
    {"name": "equal_negative_slopes", "left_slope_target": -1.0e-6, "right_slope_target": -1.0e-6},
]
NEGATIVE_EXPORT_FLAGS = [
    "apd_boundary_nullspace_sobolev_source_exported",
    "apd_nullspace_gamma_selector_source_exported",
    "apd_discrete_grid_measure_source_exported",
    "strict_dynamical_source_for_A_P_D_exported",
    "strict_phase_frequency_source_exported",
    "strict_damping_beta_eta_source_exported",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_certificate",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2576|S1526|APD nullspace Sobolev|boundary nullspace Sobolev",
        "intended_research_nonduplication": "APD.*nullspace.*roughness|nullspace.*APD.*roughness|APD.*gamma selector|gamma selector.*APD|boundary-preserving.*APD.*selector|APD.*nullspace selector dependence",
        "apd_precursors": "P2416|S1366|P2574|S1524|P2575|S1525|APD.*boundary nullspace|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def interpolation_polynomials(rows: list[dict[str, Any]]) -> tuple[Polynomial, list[Polynomial]]:
    x_nodes = np.array([row["d"] for row in rows], dtype=float)
    y_nodes = np.array([row["apd_value"] for row in rows], dtype=float)
    base = Polynomial.fit(x_nodes, y_nodes, deg=11, domain=[0.0, 11.0]).convert()
    vanish = Polynomial([1.0])
    x_poly = Polynomial([0.0, 1.0])
    for d in DOMAIN:
        vanish = vanish * Polynomial([-float(d), 1.0])
    return base, [vanish, x_poly * vanish, x_poly * x_poly * vanish]


def polynomial_combination(coefficients: np.ndarray, basis: list[Polynomial]) -> Polynomial:
    total = Polynomial([0.0])
    for coefficient, poly in zip(coefficients, basis):
        total = total + float(coefficient) * poly
    return total


def boundary_matrix(basis: list[Polynomial]) -> np.ndarray:
    return np.array([[poly.deriv(1)(0.0) for poly in basis], [poly.deriv(1)(11.0) for poly in basis]], dtype=float)


def target_solution(base: Polynomial, basis: list[Polynomial], matrix: np.ndarray, target: dict[str, Any]) -> tuple[Polynomial, np.ndarray]:
    rhs = np.array([target["left_slope_target"] - base.deriv(1)(0.0), target["right_slope_target"] - base.deriv(1)(11.0)], dtype=float)
    coeffs = np.linalg.lstsq(matrix, rhs, rcond=None)[0]
    return base + polynomial_combination(coeffs, basis), coeffs


def order_row(order: int, q0: Polynomial, null_poly: Polynomial, rows: list[dict[str, Any]], target: dict[str, Any]) -> dict[str, Any]:
    grid = np.array(DISCRETE_GRID, dtype=float)
    q_values = q0.deriv(order)(grid)
    n_values = null_poly.deriv(order)(grid)
    a_quad = float(np.dot(n_values, n_values))
    b_cross = float(np.dot(q_values, n_values))
    gamma_star = -b_cross / a_quad if a_quad > 0.0 else 0.0
    selected = q0 + gamma_star * null_poly
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    node_residuals = selected(node_points) - node_values
    left_residual = selected.deriv(1)(0.0) - target["left_slope_target"]
    right_residual = selected.deriv(1)(11.0) - target["right_slope_target"]
    return {
        "derivative_order": order,
        "objective": f"J_{order}(gamma)=sum_grid (d^{order}/dx^{order}(q0+gamma*N))^2",
        "quadratic_A": a_quad,
        "linear_cross_B": b_cross,
        "gamma_star": float(gamma_star),
        "max_abs_node_residual_at_gamma_star": float(np.max(np.abs(node_residuals))),
        "max_abs_boundary_slope_residual_at_gamma_star": float(max(abs(left_residual), abs(right_residual))),
    }


def target_row(target: dict[str, Any], base: Polynomial, basis: list[Polynomial], matrix: np.ndarray, null_vector: np.ndarray, rows: list[dict[str, Any]]) -> dict[str, Any]:
    q0, coeffs = target_solution(base, basis, matrix, target)
    null_poly = polynomial_combination(null_vector, basis)
    order_rows = [order_row(order, q0, null_poly, rows, target) for order in DERIVATIVE_ORDERS]
    rounded = {round(row["gamma_star"], 20) for row in order_rows}
    return {
        "target_name": target["name"],
        "minimum_norm_solution_coefficients_for_V_xV_x2V": [float(value) for value in coeffs],
        "order_rows": order_rows,
        "distinct_gamma_minimizers_after_rounding_1e_minus_20": len(rounded),
        "orders_audited": DERIVATIVE_ORDERS,
        "all_gamma_minimizers_preserve_apd_nodes": all(row["max_abs_node_residual_at_gamma_star"] <= 1.0e-6 for row in order_rows),
        "all_gamma_minimizers_preserve_boundary_targets": all(row["max_abs_boundary_slope_residual_at_gamma_star"] <= 1.0e-6 for row in order_rows),
        "sobolev_order_changes_nullspace_selector": len(rounded) > 1,
    }


def nullspace_sobolev_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, basis = interpolation_polynomials(rows)
    matrix = boundary_matrix(basis)
    _, singular_values, vh = np.linalg.svd(matrix)
    null_vector = vh[-1]
    target_rows = [target_row(target, base, basis, matrix, null_vector, rows) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "domain": DOMAIN,
        "discrete_grid": DISCRETE_GRID,
        "derivative_orders": DERIVATIVE_ORDERS,
        "augmented_vanishing_basis": ["V", "x*V", "x^2*V"],
        "boundary_matrix_rank": int(np.linalg.matrix_rank(matrix)),
        "boundary_matrix_nullity": int(len(basis) - np.linalg.matrix_rank(matrix)),
        "boundary_matrix_singular_values": [float(value) for value in singular_values],
        "null_vector_for_V_xV_x2V": [float(value) for value in null_vector],
        "target_rows": target_rows,
        "target_count": len(target_rows),
        "all_targets_have_order_dependent_gamma_selector": all(row["sobolev_order_changes_nullspace_selector"] for row in target_rows),
        "all_selected_gammas_preserve_nodes_and_boundaries": all(row["all_gamma_minimizers_preserve_apd_nodes"] and row["all_gamma_minimizers_preserve_boundary_targets"] for row in target_rows),
        "nullspace_sobolev_selector_remains_metric_dependent": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2574_payload = load_json(SOURCE_FILES["P2574_APD_TWO_ENDPOINT_COMPATIBILITY"])
    p2575_payload = load_json(SOURCE_FILES["P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE"])
    p2574 = theorem(p2574_payload, "apd_two_endpoint_boundary_compatibility_obstruction_audit")
    p2575 = theorem(p2575_payload, "apd_augmented_boundary_nullspace_nonuniqueness_audit")
    rows = apd_rows(p2416_payload)
    audit = nullspace_sobolev_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2576_T1_apd_boundary_nullspace_discrete_sobolev_selector_dependence_audit",
        "audited_chain": ["P2416/S1366", "P2574/S1524", "P2575/S1525"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "select the P2575 boundary-preserving nullspace gamma by minimizing discrete Sobolev roughness on a fixed grid",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2574_two_endpoint_obstruction_inherited": p2574.get("finite_apd_values_do_not_supply_compatible_two_endpoint_boundary_law") is True,
        "p2575_augmented_nullspace_inherited": p2575.get("but_augmented_boundary_solution_is_nonunique") is True,
        "apd_node_rows": rows,
        "boundary_nullspace_discrete_sobolev_selector_audit": audit,
        "finite_apd_values_and_boundary_targets_do_not_select_nullspace_gamma": audit["all_targets_have_order_dependent_gamma_selector"],
        "sobolev_nullspace_selector_is_conditional_not_strict_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not select the P2575 boundary-preserving nullspace by an unsourced roughness order or grid measure. The next honest step is to derive the APD nullspace selector, measure, and derivative order from the strict action; otherwise augmented boundary solvability remains conditional and nonunique."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2574_two_endpoint_obstruction_inherited": theorem_export["p2574_two_endpoint_obstruction_inherited"],
        "p2575_augmented_nullspace_inherited": theorem_export["p2575_augmented_nullspace_inherited"],
        "boundary_matrix_rank_two": audit["boundary_matrix_rank"] == 2,
        "boundary_matrix_nullity_one": audit["boundary_matrix_nullity"] == 1,
        "three_targets_audited": audit["target_count"] == 3,
        "all_targets_order_dependent": audit["all_targets_have_order_dependent_gamma_selector"],
        "all_selected_gammas_preserve_nodes_and_boundaries": audit["all_selected_gammas_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2576",
        "stage_id": "S1526",
        "status": "P2576_APD_BOUNDARY_NULLSPACE_DISCRETE_SOBOLEV_SELECTOR_DEPENDENCE_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_boundary_nullspace_discrete_sobolev_selector_dependence_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2574_APD_TWO_ENDPOINT_COMPATIBILITY": sha256_json(p2574_payload),
                "P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE": sha256_json(p2575_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_boundary_nullspace_discrete_sobolev_selector_dependence_audit"]["theorem_export"]
    audit = t["boundary_nullspace_discrete_sobolev_selector_audit"]
    lines = [
        "# P2576/S1526 APD boundary-nullspace discrete Sobolev selector dependence audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Derivative orders audited: `{audit['derivative_orders']}`.",
        f"- Boundary matrix rank/nullity: `{audit['boundary_matrix_rank']}/{audit['boundary_matrix_nullity']}`.",
        f"- Targets audited: `{audit['target_count']}`.",
        f"- All targets have order-dependent gamma selector: `{audit['all_targets_have_order_dependent_gamma_selector']}`.",
        f"- Selected gammas preserve nodes and boundaries: `{audit['all_selected_gammas_preserve_nodes_and_boundaries']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "On the P2575 augmented boundary family, the boundary-preserving nullspace can be parameterized by `gamma`.  Minimizing discrete Sobolev roughness over the same fixed grid chooses different `gamma` values for different derivative orders, while preserving the finite APD nodes and endpoint slopes.  Therefore the nullspace selector needs its own strict source.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD boundary-nullspace Sobolev source, gamma selector source, grid-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_boundary_nullspace_discrete_sobolev_selector_dependence_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2576/S1526` tests whether the P2575 boundary-preserving nullspace can be selected by a discrete Sobolev roughness rule.  On a fixed grid and for derivative orders `{0,1,2,3,4,12}`, the minimizing nullspace parameter `gamma` changes with the derivative order while preserving finite APD nodes and the audited endpoint slopes.  Thus even after augmented boundary solvability, the nullspace selector remains a separate strict source obligation.
""".strip()
    lag_section = """
`P2576/S1526` blocks adding a boundary-nullspace Sobolev selector to role-bearing `L_total` without sourcing its derivative order and grid/measure.  The APD nullspace can be selected in multiple metric-dependent ways; strict action data remain required.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2576/S1526 APD boundary-nullspace discrete Sobolev selector dependence guard", "## P2576/S1526 APD boundary-nullspace discrete Sobolev selector dependence guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2576/S1526 APD boundary-nullspace discrete Sobolev selector dependence Ltotal guard", "## P2576/S1526 APD boundary-nullspace discrete Sobolev selector dependence Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
