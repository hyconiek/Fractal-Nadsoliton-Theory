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
from p2578_s1528_apd_augmented_boundary_basis_metric_dependence_audit import BOUNDARY_TARGETS, DOMAIN, PROBE_POINTS, base_and_vanish, polynomial_combination
from p2580_s1530_apd_inner_product_basis_covariance_requirement_certificate import boundary_matrix, constrained_minimizer
from p2581_s1531_apd_gram_measure_moment_dependence_audit import gram_metric_for_measure, reference_basis

GEN = ROOT / "generated"
OUT = GEN / "p2582_s1532_apd_low_order_moment_matched_measure_nonuniqueness_audit.json"
MD = GEN / "p2582_s1532_apd_low_order_moment_matched_measure_nonuniqueness_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2580_APD_BASIS_COVARIANCE_REQUIREMENT": GEN / "p2580_s1530_apd_inner_product_basis_covariance_requirement_certificate.json",
    "P2581_APD_GRAM_MEASURE_DEPENDENCE": GEN / "p2581_s1531_apd_gram_measure_moment_dependence_audit.json",
}
SUPPORT_POINTS = [0.25, 3.25, 7.25, 10.75]
BASE_WEIGHTS = [1.0, 1.0, 1.0, 1.0]
LAMBDA_VALUES = [-1.0, 0.0, 1.0]
MOMENT_ORDERS_MATCHED = [0, 1, 2]
NEGATIVE_EXPORT_FLAGS = [
    "apd_low_order_moment_law_source_exported",
    "apd_moment_matched_measure_source_exported",
    "apd_positive_measure_source_exported",
    "apd_l2_inner_product_source_exported",
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
        "new_packet": "P2582|S1532|APD moment matching|moment matching.*APD",
        "intended_research_nonduplication": "APD.*moment-matched measure|moment-matched.*APD|APD.*low-order moments|low-order moments.*APD|APD.*measure continuum|measure continuum.*APD|APD.*same moments.*measure|same moments.*APD|APD.*moment nullspace|moment nullspace.*APD",
        "apd_precursors": "P2416|S1366|P2580|S1530|P2581|S1531|APD.*Gram measure|APD.*positive measure|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def moment_matrix(points: np.ndarray) -> np.ndarray:
    return np.vstack([points ** order for order in MOMENT_ORDERS_MATCHED])


def normalized_moment_null_vector(points: np.ndarray) -> np.ndarray:
    _, _, vh = np.linalg.svd(moment_matrix(points))
    vector = vh[-1]
    norm = np.linalg.norm(vector)
    if norm == 0.0:
        raise ValueError("moment matrix has no nonzero null vector")
    return vector / norm


def matched_measure_variants() -> list[dict[str, Any]]:
    points = np.array(SUPPORT_POINTS, dtype=float)
    base = np.array(BASE_WEIGHTS, dtype=float)
    null_vector = normalized_moment_null_vector(points)
    variants = []
    reference_moments = moment_matrix(points) @ base
    for lam in LAMBDA_VALUES:
        weights = base + lam * null_vector
        moments = moment_matrix(points) @ weights
        variants.append({
            "name": f"moment_matched_lambda_{lam:+.1f}".replace("+", "plus_").replace("-", "minus_").replace(".", "p"),
            "lambda": lam,
            "points": [float(value) for value in points],
            "weights": [float(value) for value in weights],
            "matched_moments_order_0_1_2": [float(value) for value in moments],
            "max_abs_moment_deviation_from_reference": float(np.max(np.abs(moments - reference_moments))),
            "all_weights_positive": bool(np.min(weights) > 0.0),
        })
    return variants


def measure_row(measure: dict[str, Any], base: Any, basis: list[Any], matrix: np.ndarray, target: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    rhs = np.array([target["left_slope_target"] - base.deriv(1)(0.0), target["right_slope_target"] - base.deriv(1)(11.0)], dtype=float)
    gram = gram_metric_for_measure(basis, measure)
    coefficients = constrained_minimizer(matrix, rhs, gram)
    selected = base + polynomial_combination(coefficients, basis)
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    node_residuals = selected(node_points) - node_values
    left_residual = selected.deriv(1)(0.0) - target["left_slope_target"]
    right_residual = selected.deriv(1)(11.0) - target["right_slope_target"]
    eigenvalues = np.linalg.eigvalsh(gram)
    return {
        "measure_name": measure["name"],
        "lambda": measure["lambda"],
        "weights": measure["weights"],
        "matched_moments_order_0_1_2": measure["matched_moments_order_0_1_2"],
        "max_abs_moment_deviation_from_reference": measure["max_abs_moment_deviation_from_reference"],
        "gram_metric_eigenvalues": [float(value) for value in eigenvalues],
        "gram_metric_is_positive_definite": bool(np.min(eigenvalues) > 0.0),
        "selected_coefficients": [float(value) for value in coefficients],
        "probe_values": [float(selected(point)) for point in PROBE_POINTS],
        "max_abs_node_residual": float(np.max(np.abs(node_residuals))),
        "max_abs_boundary_slope_residual": float(max(abs(left_residual), abs(right_residual))),
    }


def target_row(target: dict[str, Any], base: Any, basis: list[Any], matrix: np.ndarray, rows: list[dict[str, Any]], measures: list[dict[str, Any]]) -> dict[str, Any]:
    measure_rows = [measure_row(measure, base, basis, matrix, target, rows) for measure in measures]
    middle_probe_values = {round(row["probe_values"][1], 12) for row in measure_rows}
    return {
        "target_name": target["name"],
        "measure_rows": measure_rows,
        "moment_matched_measure_count": len(measure_rows),
        "distinct_middle_probe_values_after_rounding_1e_minus_12": len(middle_probe_values),
        "all_measures_share_low_order_moments": all(row["max_abs_moment_deviation_from_reference"] <= 1.0e-12 for row in measure_rows),
        "all_measures_positive": all(min(row["weights"]) > 0.0 for row in measure_rows),
        "all_gram_metrics_positive_definite": all(row["gram_metric_is_positive_definite"] for row in measure_rows),
        "all_minimizers_preserve_apd_nodes": all(row["max_abs_node_residual"] <= 1.0e-6 for row in measure_rows),
        "all_minimizers_satisfy_boundary_targets": all(row["max_abs_boundary_slope_residual"] <= 1.0e-6 for row in measure_rows),
        "low_order_moments_do_not_select_apd_measure": len(middle_probe_values) > 1,
    }


def moment_matching_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    basis = reference_basis(vanish)
    matrix = boundary_matrix(basis)
    measures = matched_measure_variants()
    points = np.array(SUPPORT_POINTS, dtype=float)
    null_vector = normalized_moment_null_vector(points)
    target_rows = [target_row(target, base, basis, matrix, rows, measures) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "support_points": SUPPORT_POINTS,
        "base_weights": BASE_WEIGHTS,
        "moment_orders_matched": MOMENT_ORDERS_MATCHED,
        "moment_null_vector": [float(value) for value in null_vector],
        "lambda_values": LAMBDA_VALUES,
        "matched_measure_variants": measures,
        "target_rows": target_rows,
        "target_count": len(target_rows),
        "moment_matched_measure_count": len(measures),
        "all_measures_positive": all(measure["all_weights_positive"] for measure in measures),
        "all_measures_share_low_order_moments": all(measure["max_abs_moment_deviation_from_reference"] <= 1.0e-12 for measure in measures),
        "all_targets_low_order_moment_nonunique": all(row["low_order_moments_do_not_select_apd_measure"] for row in target_rows),
        "all_gram_metrics_positive_definite": all(row["all_gram_metrics_positive_definite"] for row in target_rows),
        "all_minimizers_preserve_nodes_and_boundaries": all(row["all_minimizers_preserve_apd_nodes"] and row["all_minimizers_satisfy_boundary_targets"] for row in target_rows),
        "low_order_moment_law_is_unsourced_selector": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2580_payload = load_json(SOURCE_FILES["P2580_APD_BASIS_COVARIANCE_REQUIREMENT"])
    p2581_payload = load_json(SOURCE_FILES["P2581_APD_GRAM_MEASURE_DEPENDENCE"])
    p2580 = theorem(p2580_payload, "apd_inner_product_basis_covariance_requirement_certificate")
    p2581 = theorem(p2581_payload, "apd_gram_measure_moment_dependence_audit")
    rows = apd_rows(p2416_payload)
    audit = moment_matching_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2582_T1_apd_low_order_moment_matched_measure_nonuniqueness_audit",
        "audited_chain": ["P2416/S1366", "P2580/S1530", "P2581/S1531"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "select the APD Gram measure by fixing low-order measure moments mass/mean/second raw moment",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2580_metric_covariance_requirement_inherited": p2580.get("coordinate_covariance_removes_basis_artifact_but_not_metric_choice") is True,
        "p2581_measure_dependence_inherited": p2581.get("finite_apd_values_and_boundary_targets_do_not_select_positive_measure") is True,
        "apd_node_rows": rows,
        "apd_low_order_moment_matched_measure_nonuniqueness_audit": audit,
        "finite_low_order_moments_do_not_select_positive_measure": audit["all_targets_low_order_moment_nonunique"],
        "low_order_moment_selector_is_conditional_not_strict_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not promote low-order APD measure moments into a strict source. P2582 constructs positive measures with identical mass, first moment, and second raw moment that still induce different Gram metrics and off-node APD dynamics. The next honest step is to derive the full APD measure density or all moment constraints from strict nadsoliton dynamics, rather than fixing only a finite moment prefix."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2580_metric_covariance_requirement_inherited": theorem_export["p2580_metric_covariance_requirement_inherited"],
        "p2581_measure_dependence_inherited": theorem_export["p2581_measure_dependence_inherited"],
        "three_targets_audited": audit["target_count"] == 3,
        "three_moment_matched_measures": audit["moment_matched_measure_count"] == 3,
        "all_measures_positive": audit["all_measures_positive"],
        "all_measures_share_low_order_moments": audit["all_measures_share_low_order_moments"],
        "all_targets_low_order_moment_nonunique": audit["all_targets_low_order_moment_nonunique"],
        "all_gram_metrics_positive_definite": audit["all_gram_metrics_positive_definite"],
        "all_minimizers_preserve_nodes_and_boundaries": audit["all_minimizers_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2582",
        "stage_id": "S1532",
        "status": "P2582_APD_LOW_ORDER_MOMENT_MATCHED_MEASURE_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_low_order_moment_matched_measure_nonuniqueness_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2580_APD_BASIS_COVARIANCE_REQUIREMENT": sha256_json(p2580_payload),
                "P2581_APD_GRAM_MEASURE_DEPENDENCE": sha256_json(p2581_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_low_order_moment_matched_measure_nonuniqueness_audit"]["theorem_export"]
    audit = t["apd_low_order_moment_matched_measure_nonuniqueness_audit"]
    lines = [
        "# P2582/S1532 APD low-order moment-matched measure nonuniqueness audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Boundary targets audited: `{audit['target_count']}`.",
        f"- Moment-matched positive measures audited: `{audit['moment_matched_measure_count']}`.",
        f"- Moment orders matched: `{audit['moment_orders_matched']}`.",
        f"- All targets low-order moment nonunique: `{audit['all_targets_low_order_moment_nonunique']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "Matching mass, first moment, and second raw moment does not determine the APD Gram measure.  P2582 constructs a positive moment-nullspace family with identical low-order moments; every audited measure preserves the APD nodes and endpoint slopes, but the selected off-node APD dynamics still changes.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No low-order APD moment-law source, moment-matched measure source, positive-measure source, L2-inner-product source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_low_order_moment_matched_measure_nonuniqueness_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2582/S1532 APD low-order moment-matched measure nonuniqueness guard

`P2582/S1532` tests whether the P2581 measure dependence can be removed by fixing low-order moments of the positive APD measure.  A three-measure positive family on the same support has identical mass, first moment, and second raw moment, and each induced Gram metric preserves the finite APD nodes plus endpoint slopes; nevertheless the selected off-node APD dynamics still varies.  Thus a finite low-order moment law is not `strict_dynamical_source_for_A_P_D`.
""".strip()
    lag_section = """
## P2582/S1532 APD low-order moment-matched measure nonuniqueness Ltotal guard

`P2582/S1532` blocks adding a role-bearing APD Gram kinetic term to `L_total` by specifying only a finite prefix of measure moments.  The strict action must source the full measure density or a complete moment law; matching mass/mean/variance-like data still leaves measure-dependent APD dynamics.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2582/S1532 APD low-order moment-matched measure nonuniqueness guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2582/S1532 APD low-order moment-matched measure nonuniqueness Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
