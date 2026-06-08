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
OUT = GEN / "p2575_s1525_apd_augmented_boundary_nullspace_nonuniqueness_audit.json"
MD = GEN / "p2575_s1525_apd_augmented_boundary_nullspace_nonuniqueness_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2573_APD_INVERSE_BOUNDARY_TUNABILITY": GEN / "p2573_s1523_apd_boundary_penalty_inverse_target_tunability_audit.json",
    "P2574_APD_TWO_ENDPOINT_COMPATIBILITY": GEN / "p2574_s1524_apd_two_endpoint_boundary_compatibility_obstruction_audit.json",
}
DOMAIN = list(range(12))
PROBE_POINTS = [0.5, 5.5, 11.5]
NULLSPACE_GAMMA_VALUES = [0.0, 1.0e-14, -1.0e-14]
BOUNDARY_TARGETS = [
    {"name": "zero_zero_neumann", "left_slope_target": 0.0, "right_slope_target": 0.0},
    {"name": "equal_positive_slopes", "left_slope_target": 1.0e-6, "right_slope_target": 1.0e-6},
    {"name": "equal_negative_slopes", "left_slope_target": -1.0e-6, "right_slope_target": -1.0e-6},
]
NEGATIVE_EXPORT_FLAGS = [
    "apd_augmented_boundary_source_exported",
    "apd_vanishing_mode_source_exported",
    "apd_boundary_nullspace_selector_source_exported",
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
        "new_packet": "P2575|S1525|APD augmented boundary|two-endpoint augmented",
        "intended_research_nonduplication": "APD.*augmented.*boundary|augmented.*APD.*boundary|APD.*boundary nullspace|boundary nullspace.*APD|APD.*multi-vanishing|multi-vanishing.*APD|two-endpoint.*nonuniqueness.*APD|APD.*boundary.*nonuniqueness",
        "apd_precursors": "P2416|S1366|P2573|S1523|P2574|S1524|APD.*boundary|APD.*Neumann|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def interpolation_polynomials(rows: list[dict[str, Any]]) -> tuple[Polynomial, Polynomial, list[Polynomial]]:
    x_nodes = np.array([row["d"] for row in rows], dtype=float)
    y_nodes = np.array([row["apd_value"] for row in rows], dtype=float)
    base = Polynomial.fit(x_nodes, y_nodes, deg=11, domain=[0.0, 11.0]).convert()
    vanish = Polynomial([1.0])
    x_poly = Polynomial([0.0, 1.0])
    for d in DOMAIN:
        vanish = vanish * Polynomial([-float(d), 1.0])
    basis = [vanish, x_poly * vanish, x_poly * x_poly * vanish]
    return base, vanish, basis


def polynomial_combination(coefficients: np.ndarray, basis: list[Polynomial]) -> Polynomial:
    total = Polynomial([0.0])
    for coefficient, poly in zip(coefficients, basis):
        total = total + float(coefficient) * poly
    return total


def boundary_matrix(basis: list[Polynomial]) -> np.ndarray:
    return np.array([[poly.deriv(1)(0.0) for poly in basis], [poly.deriv(1)(11.0) for poly in basis]], dtype=float)


def target_row(target: dict[str, Any], base: Polynomial, basis: list[Polynomial], matrix: np.ndarray, null_vector: np.ndarray, rows: list[dict[str, Any]]) -> dict[str, Any]:
    rhs = np.array([target["left_slope_target"] - base.deriv(1)(0.0), target["right_slope_target"] - base.deriv(1)(11.0)], dtype=float)
    minimum_norm_solution = np.linalg.lstsq(matrix, rhs, rcond=None)[0]
    family_rows = []
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    for gamma in NULLSPACE_GAMMA_VALUES:
        coefficients = minimum_norm_solution + gamma * null_vector
        correction = polynomial_combination(coefficients, basis)
        selected = base + correction
        node_residuals = selected(node_points) - node_values
        left_residual = selected.deriv(1)(0.0) - target["left_slope_target"]
        right_residual = selected.deriv(1)(11.0) - target["right_slope_target"]
        probe_values = [float(selected(point)) for point in PROBE_POINTS]
        base_probe_values = [float(base(point)) for point in PROBE_POINTS]
        family_rows.append({
            "gamma": gamma,
            "coefficients_for_V_xV_x2V": [float(value) for value in coefficients],
            "max_abs_node_residual": float(np.max(np.abs(node_residuals))),
            "max_abs_boundary_slope_residual": float(max(abs(left_residual), abs(right_residual))),
            "probe_values": probe_values,
            "probe_deltas_from_base": [float(value - base_value) for value, base_value in zip(probe_values, base_probe_values)],
        })
    nonzero = [row for row in family_rows if row["gamma"] != 0.0]
    return {
        "target_name": target["name"],
        "left_slope_target": target["left_slope_target"],
        "right_slope_target": target["right_slope_target"],
        "minimum_norm_solution_coefficients_for_V_xV_x2V": [float(value) for value in minimum_norm_solution],
        "family_rows": family_rows,
        "all_family_rows_preserve_apd_nodes": all(row["max_abs_node_residual"] <= 1.0e-6 for row in family_rows),
        "all_family_rows_satisfy_boundary_targets": all(row["max_abs_boundary_slope_residual"] <= 1.0e-8 for row in family_rows),
        "nonzero_nullspace_members_change_off_node_values": all(any(abs(value) > 0.0 for value in row["probe_deltas_from_base"]) for row in nonzero),
    }


def augmented_boundary_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish, basis = interpolation_polynomials(rows)
    matrix = boundary_matrix(basis)
    _, singular_values, vh = np.linalg.svd(matrix)
    null_vector = vh[-1]
    target_rows = [target_row(target, base, basis, matrix, null_vector, rows) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "domain": DOMAIN,
        "probe_points": PROBE_POINTS,
        "augmented_vanishing_basis": ["V=prod_{d=0}^{11}(x-d)", "x*V", "x^2*V"],
        "boundary_matrix_rows_left_right_columns_V_xV_x2V": matrix.tolist(),
        "boundary_matrix_rank": int(np.linalg.matrix_rank(matrix)),
        "boundary_matrix_nullity": int(len(basis) - np.linalg.matrix_rank(matrix)),
        "boundary_matrix_singular_values": [float(value) for value in singular_values],
        "null_vector_for_V_xV_x2V": [float(value) for value in null_vector],
        "target_rows": target_rows,
        "target_count": len(target_rows),
        "all_targets_solved_by_augmented_basis": all(row["all_family_rows_satisfy_boundary_targets"] for row in target_rows),
        "all_targets_preserve_finite_apd_nodes": all(row["all_family_rows_preserve_apd_nodes"] for row in target_rows),
        "nullspace_keeps_boundary_targets_but_changes_off_node_values": all(row["nonzero_nullspace_members_change_off_node_values"] for row in target_rows),
        "augmented_boundary_solution_is_not_unique": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2573_payload = load_json(SOURCE_FILES["P2573_APD_INVERSE_BOUNDARY_TUNABILITY"])
    p2574_payload = load_json(SOURCE_FILES["P2574_APD_TWO_ENDPOINT_COMPATIBILITY"])
    p2573 = theorem(p2573_payload, "apd_boundary_penalty_inverse_target_tunability_audit")
    p2574 = theorem(p2574_payload, "apd_two_endpoint_boundary_compatibility_obstruction_audit")
    rows = apd_rows(p2416_payload)
    audit = augmented_boundary_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2575_T1_apd_augmented_boundary_nullspace_nonuniqueness_audit",
        "audited_chain": ["P2416/S1366", "P2573/S1523", "P2574/S1524"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "augment the P2574 one-parameter family with x*V and x^2*V to satisfy two endpoint boundary targets",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2573_inverse_boundary_tunability_inherited": p2573.get("finite_apd_values_do_not_select_target_or_penalty_law") is True,
        "p2574_two_endpoint_obstruction_inherited": p2574.get("finite_apd_values_do_not_supply_compatible_two_endpoint_boundary_law") is True,
        "apd_node_rows": rows,
        "augmented_boundary_nullspace_audit": audit,
        "two_endpoint_boundary_can_be_satisfied_after_adding_vanishing_modes": audit["all_targets_solved_by_augmented_basis"],
        "but_augmented_boundary_solution_is_nonunique": audit["augmented_boundary_solution_is_not_unique"],
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not repair two-endpoint APD incompatibility by freely adding extra vanishing modes and then claim a source theorem. The next honest step is to derive the admissible APD function space and boundary data from the strict action; otherwise the boundary problem can be made solvable but remains nonunique through a nullspace."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2573_inverse_boundary_tunability_inherited": theorem_export["p2573_inverse_boundary_tunability_inherited"],
        "p2574_two_endpoint_obstruction_inherited": theorem_export["p2574_two_endpoint_obstruction_inherited"],
        "twelve_apd_rows": len(rows) == 12,
        "boundary_matrix_rank_two": audit["boundary_matrix_rank"] == 2,
        "boundary_matrix_nullity_one": audit["boundary_matrix_nullity"] == 1,
        "all_targets_solved_by_augmented_basis": audit["all_targets_solved_by_augmented_basis"],
        "nullspace_changes_off_node_values": audit["nullspace_keeps_boundary_targets_but_changes_off_node_values"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2575",
        "stage_id": "S1525",
        "status": "P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_augmented_boundary_nullspace_nonuniqueness_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2573_APD_INVERSE_BOUNDARY_TUNABILITY": sha256_json(p2573_payload),
                "P2574_APD_TWO_ENDPOINT_COMPATIBILITY": sha256_json(p2574_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_augmented_boundary_nullspace_nonuniqueness_audit"]["theorem_export"]
    audit = t["augmented_boundary_nullspace_audit"]
    lines = [
        "# P2575/S1525 APD augmented-boundary nullspace nonuniqueness audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Augmented basis: `{audit['augmented_vanishing_basis']}`.",
        f"- Boundary matrix rank/nullity: `{audit['boundary_matrix_rank']}/{audit['boundary_matrix_nullity']}`.",
        f"- Targets audited: `{audit['target_count']}`.",
        f"- All targets solved by augmented basis: `{audit['all_targets_solved_by_augmented_basis']}`.",
        f"- Nullspace changes off-node values while preserving boundaries: `{audit['nullspace_keeps_boundary_targets_but_changes_off_node_values']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "Adding `V`, `x*V`, and `x^2*V` repairs the P2574 two-endpoint compatibility obstruction, but the resulting boundary matrix has rank `2` and nullity `1`.  Therefore the augmented family can satisfy the audited endpoint slopes while still carrying a null direction that preserves all finite APD nodes and boundary slopes but changes off-node values.  Solvability by extra vanishing modes is not a source theorem.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD augmented boundary source, vanishing-mode source, boundary nullspace selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_augmented_boundary_nullspace_nonuniqueness_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2575/S1525` tests the repair move after P2574: add more node-vanishing modes `V`, `x*V`, `x^2*V` so two endpoint APD slope targets become solvable.  The boundary matrix has rank `2` and nullity `1`; all audited targets can be met while preserving finite APD nodes, but the null direction still changes off-node values.  Thus augmented boundary solvability is not `strict_dynamical_source_for_A_P_D` without a sourced admissible function space.
""".strip()
    lag_section = """
`P2575/S1525` blocks repairing APD boundary incompatibility inside role-bearing `L_total` by freely adding vanishing modes.  Extra modes can satisfy endpoint data but leave a boundary-preserving nullspace; the strict action must source the allowed APD function space and nullspace selector.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2575/S1525 APD augmented-boundary nullspace nonuniqueness guard", "## P2575/S1525 APD augmented-boundary nullspace nonuniqueness guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2575/S1525 APD augmented-boundary nullspace nonuniqueness Ltotal guard", "## P2575/S1525 APD augmented-boundary nullspace nonuniqueness Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
