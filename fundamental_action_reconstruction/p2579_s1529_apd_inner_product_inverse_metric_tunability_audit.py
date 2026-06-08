#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import numpy as np

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2569_s1519_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate import apd_rows
from p2575_s1525_apd_augmented_boundary_nullspace_nonuniqueness_audit import BOUNDARY_TARGETS, DOMAIN, PROBE_POINTS, boundary_matrix, interpolation_polynomials, polynomial_combination

GEN = ROOT / "generated"
OUT = GEN / "p2579_s1529_apd_inner_product_inverse_metric_tunability_audit.json"
MD = GEN / "p2579_s1529_apd_inner_product_inverse_metric_tunability_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE": GEN / "p2575_s1525_apd_augmented_boundary_nullspace_nonuniqueness_audit.json",
    "P2578_APD_BASIS_METRIC_DEPENDENCE": GEN / "p2578_s1528_apd_augmented_boundary_basis_metric_dependence_audit.json",
}
TARGET_GAMMA_VALUES = [-2.0e-14, -1.0e-14, 0.0, 1.0e-14, 2.0e-14]
NEGATIVE_EXPORT_FLAGS = [
    "apd_function_space_inner_product_source_exported",
    "apd_spd_metric_source_exported",
    "apd_inverse_metric_selector_source_exported",
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
        "new_packet": "P2579|S1529|APD inner product|inner product.*APD",
        "intended_research_nonduplication": "APD.*metric tunability|metric tunability.*APD|SPD metric.*APD|APD.*SPD metric|APD.*quadratic metric|quadratic metric.*APD|APD.*inverse metric|inverse metric.*APD|function-space inner product.*APD|APD.*target gamma.*metric",
        "apd_precursors": "P2416|S1366|P2575|S1525|P2578|S1528|APD.*basis metric|APD.*nullspace|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def null_vector_from_matrix(matrix: np.ndarray) -> np.ndarray:
    _, _, vh = np.linalg.svd(matrix)
    vector = vh[-1]
    norm = np.linalg.norm(vector)
    if norm == 0.0:
        raise ValueError("boundary matrix has no computable null vector")
    return vector / norm


def spd_metric_for_target_gamma(minimum_solution: np.ndarray, null_vector: np.ndarray, target_gamma: float) -> np.ndarray:
    row_direction_norm = np.linalg.norm(minimum_solution)
    if row_direction_norm == 0.0:
        raise ValueError("minimum solution is zero; inverse metric witness would be degenerate")
    n_hat = null_vector / np.linalg.norm(null_vector)
    row_hat = minimum_solution / row_direction_norm
    third = np.cross(n_hat, row_hat)
    third_norm = np.linalg.norm(third)
    if third_norm == 0.0:
        raise ValueError("minimum solution and null vector are linearly dependent")
    third = third / third_norm
    q = np.column_stack([n_hat, row_hat, third])
    off_diagonal = -target_gamma / row_direction_norm
    local_metric = np.array([
        [1.0, off_diagonal, 0.0],
        [off_diagonal, off_diagonal * off_diagonal + 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ], dtype=float)
    return q @ local_metric @ q.T


def metric_witness_row(target: dict[str, Any], base: Any, basis: list[Any], matrix: np.ndarray, null_vector: np.ndarray, rows: list[dict[str, Any]], target_gamma: float) -> dict[str, Any]:
    rhs = np.array([target["left_slope_target"] - base.deriv(1)(0.0), target["right_slope_target"] - base.deriv(1)(11.0)], dtype=float)
    minimum_solution = np.linalg.lstsq(matrix, rhs, rcond=None)[0]
    metric = spd_metric_for_target_gamma(minimum_solution, null_vector, target_gamma)
    selected_coefficients = minimum_solution + target_gamma * null_vector
    selected = base + polynomial_combination(selected_coefficients, basis)
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    node_residuals = selected(node_points) - node_values
    left_residual = selected.deriv(1)(0.0) - target["left_slope_target"]
    right_residual = selected.deriv(1)(11.0) - target["right_slope_target"]
    stationarity_residual = float(null_vector @ metric @ selected_coefficients)
    eigenvalues = np.linalg.eigvalsh(metric)
    recovered_gamma = float(np.dot(selected_coefficients - minimum_solution, null_vector) / np.dot(null_vector, null_vector))
    return {
        "target_gamma": target_gamma,
        "recovered_gamma_from_coefficients": recovered_gamma,
        "spd_metric_3x3": [[float(value) for value in row] for row in metric],
        "metric_eigenvalues": [float(value) for value in eigenvalues],
        "metric_is_positive_definite": bool(np.min(eigenvalues) > 0.0),
        "selected_coefficients": [float(value) for value in selected_coefficients],
        "probe_values": [float(selected(point)) for point in PROBE_POINTS],
        "metric_stationarity_residual": abs(stationarity_residual),
        "max_abs_node_residual": float(np.max(np.abs(node_residuals))),
        "max_abs_boundary_slope_residual": float(max(abs(left_residual), abs(right_residual))),
    }


def target_row(target: dict[str, Any], base: Any, basis: list[Any], matrix: np.ndarray, null_vector: np.ndarray, rows: list[dict[str, Any]]) -> dict[str, Any]:
    witness_rows = [metric_witness_row(target, base, basis, matrix, null_vector, rows, gamma) for gamma in TARGET_GAMMA_VALUES]
    middle_probe_values = {round(row["probe_values"][1], 18) for row in witness_rows}
    return {
        "target_name": target["name"],
        "metric_witness_rows": witness_rows,
        "target_gamma_values": TARGET_GAMMA_VALUES,
        "metric_witness_count": len(witness_rows),
        "distinct_middle_probe_values_after_rounding_1e_minus_18": len(middle_probe_values),
        "all_witness_metrics_positive_definite": all(row["metric_is_positive_definite"] for row in witness_rows),
        "all_witnesses_hit_requested_gamma": all(abs(row["recovered_gamma_from_coefficients"] - row["target_gamma"]) <= 1.0e-24 for row in witness_rows),
        "all_witnesses_stationary_for_their_metric": all(row["metric_stationarity_residual"] <= 1.0e-20 for row in witness_rows),
        "all_witnesses_preserve_apd_nodes": all(row["max_abs_node_residual"] <= 1.0e-6 for row in witness_rows),
        "all_witnesses_satisfy_boundary_targets": all(row["max_abs_boundary_slope_residual"] <= 1.0e-6 for row in witness_rows),
        "spd_inner_product_can_select_multiple_gammas": len(middle_probe_values) > 1,
    }


def inverse_metric_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, _, basis = interpolation_polynomials(rows)
    matrix = boundary_matrix(basis)
    null_vector = null_vector_from_matrix(matrix)
    target_rows = [target_row(target, base, basis, matrix, null_vector, rows) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "domain": DOMAIN,
        "probe_points": PROBE_POINTS,
        "target_gamma_values": TARGET_GAMMA_VALUES,
        "coefficient_basis": ["V", "xV", "x^2V"],
        "boundary_matrix_rank": int(np.linalg.matrix_rank(matrix)),
        "coefficient_nullity": int(matrix.shape[1] - np.linalg.matrix_rank(matrix)),
        "normalized_null_vector": [float(value) for value in null_vector],
        "target_rows": target_rows,
        "target_count": len(target_rows),
        "all_targets_inverse_metric_tunable": all(row["spd_inner_product_can_select_multiple_gammas"] for row in target_rows),
        "all_constructed_metrics_positive_definite": all(row["all_witness_metrics_positive_definite"] for row in target_rows),
        "all_constructed_metrics_preserve_nodes_boundaries_and_stationarity": all(
            row["all_witnesses_hit_requested_gamma"]
            and row["all_witnesses_stationary_for_their_metric"]
            and row["all_witnesses_preserve_apd_nodes"]
            and row["all_witnesses_satisfy_boundary_targets"]
            for row in target_rows
        ),
        "function_space_inner_product_is_unsourced_selector": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2575_payload = load_json(SOURCE_FILES["P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE"])
    p2578_payload = load_json(SOURCE_FILES["P2578_APD_BASIS_METRIC_DEPENDENCE"])
    p2575 = theorem(p2575_payload, "apd_augmented_boundary_nullspace_nonuniqueness_audit")
    p2578 = theorem(p2578_payload, "apd_augmented_boundary_basis_metric_dependence_audit")
    rows = apd_rows(p2416_payload)
    audit = inverse_metric_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2579_T1_apd_inner_product_inverse_metric_tunability_audit",
        "audited_chain": ["P2416/S1366", "P2575/S1525", "P2578/S1528"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "source-free SPD inner product minimization on the APD augmented-boundary coefficient nullspace",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2575_augmented_nullspace_inherited": p2575.get("but_augmented_boundary_solution_is_nonunique") is True,
        "p2578_basis_metric_dependence_inherited": p2578.get("finite_apd_values_and_boundary_targets_do_not_select_coordinate_metric") is True,
        "apd_node_rows": rows,
        "apd_inner_product_inverse_metric_tunability_audit": audit,
        "finite_apd_values_and_boundary_targets_do_not_select_function_space_inner_product": audit["all_targets_inverse_metric_tunable"],
        "spd_metric_minimizer_selector_is_conditional_not_strict_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not promote an APD inner-product minimizer unless the strict action derives the inner product itself. P2579 shows the inverse problem is tunable: for the same APD nodes and boundary targets, positive-definite metrics can be constructed to select different nullspace gammas. The next honest step is to derive a canonical APD kinetic/measure term from strict nadsoliton dynamics, not to choose a convenient SPD metric post hoc."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2575_augmented_nullspace_inherited": theorem_export["p2575_augmented_nullspace_inherited"],
        "p2578_basis_metric_dependence_inherited": theorem_export["p2578_basis_metric_dependence_inherited"],
        "three_targets_audited": audit["target_count"] == 3,
        "five_gamma_targets_per_boundary_target": len(audit["target_gamma_values"]) == 5,
        "rank_two_nullity_one_boundary_matrix": audit["boundary_matrix_rank"] == 2 and audit["coefficient_nullity"] == 1,
        "all_targets_inverse_metric_tunable": audit["all_targets_inverse_metric_tunable"],
        "all_constructed_metrics_positive_definite": audit["all_constructed_metrics_positive_definite"],
        "all_constructed_metrics_preserve_nodes_boundaries_and_stationarity": audit["all_constructed_metrics_preserve_nodes_boundaries_and_stationarity"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2579",
        "stage_id": "S1529",
        "status": "P2579_APD_INNER_PRODUCT_INVERSE_METRIC_TUNABILITY_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_inner_product_inverse_metric_tunability_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE": sha256_json(p2575_payload),
                "P2578_APD_BASIS_METRIC_DEPENDENCE": sha256_json(p2578_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_inner_product_inverse_metric_tunability_audit"]["theorem_export"]
    audit = t["apd_inner_product_inverse_metric_tunability_audit"]
    lines = [
        "# P2579/S1529 APD inner-product inverse metric tunability audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Boundary targets audited: `{audit['target_count']}`.",
        f"- Target gammas per boundary target: `{len(audit['target_gamma_values'])}`.",
        f"- Boundary matrix rank/nullity: `{audit['boundary_matrix_rank']}/{audit['coefficient_nullity']}`.",
        f"- All targets inverse-metric tunable: `{audit['all_targets_inverse_metric_tunable']}`.",
        f"- Constructed metrics positive definite: `{audit['all_constructed_metrics_positive_definite']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "For an affine one-dimensional APD boundary-nullspace family, the minimizer of `c^T G c` depends on the positive-definite metric `G`.  P2579 constructs SPD metrics that make prescribed nullspace gammas stationary while preserving finite APD nodes and endpoint slopes.  Thus an inner-product minimizer is only as strict as the sourced inner product.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD function-space inner-product source, SPD-metric source, inverse-metric selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_inner_product_inverse_metric_tunability_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2579/S1529 APD inner-product inverse metric tunability guard

`P2579/S1529` sharpens the P2578 basis-metric obstruction by solving the inverse APD metric problem.  In the rank-2/nullity-1 augmented-boundary coefficient family, positive-definite metrics can be constructed so that different nullspace `gamma` values are the minimizers of `c^T G c`, while finite APD nodes and endpoint slopes remain fixed.  Thus even a coordinate-free-looking SPD inner-product minimizer is not `strict_dynamical_source_for_A_P_D` unless the strict action sources the inner product itself.
""".strip()
    lag_section = """
## P2579/S1529 APD inner-product inverse metric tunability Ltotal guard

`P2579/S1529` blocks adding an APD inner-product minimizer to role-bearing `L_total` by metric choice.  Since SPD metrics can be inverse-tuned to select different APD boundary-nullspace dynamics under the same value and boundary data, the APD kinetic/measure term must be derived from strict nadsoliton dynamics before it can carry role-bearing force.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2579/S1529 APD inner-product inverse metric tunability guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2579/S1529 APD inner-product inverse metric tunability Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
