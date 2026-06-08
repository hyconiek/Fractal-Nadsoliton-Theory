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
OUT = GEN / "p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit.json"
MD = GEN / "p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2582_LOW_ORDER_MOMENT_MATCHING": GEN / "p2582_s1532_apd_low_order_moment_matched_measure_nonuniqueness_audit.json",
    "P2583_FINITE_MOMENT_PREFIX_LADDER": GEN / "p2583_s1533_apd_finite_moment_prefix_measure_ladder_audit.json",
}
SUPPORT_VARIANTS = [
    {"name": "left_center_right_four_point_support", "points": [0.25, 3.25, 7.25, 10.75], "weights": [1.0, 1.0, 1.0, 1.0]},
    {"name": "near_uniform_four_point_support", "points": [0.5, 3.5, 6.5, 10.5], "weights": [1.0, 1.0, 1.0, 1.0]},
    {"name": "edge_weighted_four_point_support", "points": [0.75, 2.75, 8.25, 10.25], "weights": [1.0, 1.0, 1.0, 1.0]},
]
FULL_MOMENT_ORDERS = [0, 1, 2, 3]
NEGATIVE_EXPORT_FLAGS = [
    "apd_finite_support_source_exported",
    "apd_full_moment_law_source_exported",
    "apd_support_selection_source_exported",
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
        "new_packet": "P2584|S1534|APD full moment|full moment.*APD",
        "intended_research_nonduplication": "APD.*finite support.*moment|finite support.*APD.*moment|APD.*support dependence|support dependence.*APD|APD.*Vandermonde moment|Vandermonde moment.*APD|APD.*conditional moment uniqueness|conditional moment uniqueness.*APD|APD.*support.*measure source|support.*measure source.*APD",
        "apd_precursors": "P2416|S1366|P2582|S1532|P2583|S1533|APD.*moment prefix|APD.*moment-matched|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def vandermonde_moment_matrix(points: np.ndarray) -> np.ndarray:
    return np.vstack([points ** order for order in FULL_MOMENT_ORDERS])


def full_moments(points: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return vandermonde_moment_matrix(points) @ weights


def support_witness(support: dict[str, Any]) -> dict[str, Any]:
    points = np.array(support["points"], dtype=float)
    weights = np.array(support["weights"], dtype=float)
    matrix = vandermonde_moment_matrix(points)
    moments = full_moments(points, weights)
    recovered_weights = np.linalg.solve(matrix, moments)
    return {
        "support_name": support["name"],
        "points": [float(value) for value in points],
        "input_weights": [float(value) for value in weights],
        "full_moment_orders": FULL_MOMENT_ORDERS,
        "full_moments": [float(value) for value in moments],
        "vandermonde_determinant": float(np.linalg.det(matrix)),
        "vandermonde_rank": int(np.linalg.matrix_rank(matrix)),
        "recovered_weights_from_full_moments": [float(value) for value in recovered_weights],
        "max_abs_recovered_weight_error": float(np.max(np.abs(recovered_weights - weights))),
        "full_moments_conditionally_select_weights_on_fixed_support": bool(np.max(np.abs(recovered_weights - weights)) <= 1.0e-10),
        "all_recovered_weights_positive": bool(np.min(recovered_weights) > 0.0),
    }


def support_measure_from_witness(witness: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": witness["support_name"],
        "points": witness["points"],
        "weights": witness["recovered_weights_from_full_moments"],
    }


def target_support_row(target: dict[str, Any], base: Any, basis: list[Any], matrix: np.ndarray, rows: list[dict[str, Any]], witness: dict[str, Any]) -> dict[str, Any]:
    measure = support_measure_from_witness(witness)
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
        "support_name": witness["support_name"],
        "gram_metric_eigenvalues": [float(value) for value in eigenvalues],
        "gram_metric_is_positive_definite": bool(np.min(eigenvalues) > 0.0),
        "selected_coefficients": [float(value) for value in coefficients],
        "probe_values": [float(selected(point)) for point in PROBE_POINTS],
        "max_abs_node_residual": float(np.max(np.abs(node_residuals))),
        "max_abs_boundary_slope_residual": float(max(abs(left_residual), abs(right_residual))),
    }


def target_row(target: dict[str, Any], base: Any, basis: list[Any], matrix: np.ndarray, rows: list[dict[str, Any]], support_witnesses: list[dict[str, Any]]) -> dict[str, Any]:
    support_rows = [target_support_row(target, base, basis, matrix, rows, witness) for witness in support_witnesses]
    middle_probe_values = {round(row["probe_values"][1], 12) for row in support_rows}
    return {
        "target_name": target["name"],
        "support_rows": support_rows,
        "support_count": len(support_rows),
        "distinct_middle_probe_values_after_rounding_1e_minus_12": len(middle_probe_values),
        "all_gram_metrics_positive_definite": all(row["gram_metric_is_positive_definite"] for row in support_rows),
        "all_minimizers_preserve_apd_nodes": all(row["max_abs_node_residual"] <= 1.0e-6 for row in support_rows),
        "all_minimizers_satisfy_boundary_targets": all(row["max_abs_boundary_slope_residual"] <= 1.0e-6 for row in support_rows),
        "fixed_support_full_moments_do_not_select_support": len(middle_probe_values) > 1,
    }


def full_moment_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    basis = reference_basis(vanish)
    matrix = boundary_matrix(basis)
    support_witnesses = [support_witness(support) for support in SUPPORT_VARIANTS]
    target_rows = [target_row(target, base, basis, matrix, rows, support_witnesses) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "full_moment_orders": FULL_MOMENT_ORDERS,
        "support_witnesses": support_witnesses,
        "target_rows": target_rows,
        "support_count": len(support_witnesses),
        "target_count": len(target_rows),
        "all_fixed_supports_have_vandermonde_rank_four": all(witness["vandermonde_rank"] == 4 for witness in support_witnesses),
        "all_fixed_support_weights_conditionally_unique": all(witness["full_moments_conditionally_select_weights_on_fixed_support"] for witness in support_witnesses),
        "all_recovered_weights_positive": all(witness["all_recovered_weights_positive"] for witness in support_witnesses),
        "all_targets_support_dependent": all(row["fixed_support_full_moments_do_not_select_support"] for row in target_rows),
        "all_gram_metrics_positive_definite": all(row["all_gram_metrics_positive_definite"] for row in target_rows),
        "all_minimizers_preserve_nodes_and_boundaries": all(row["all_minimizers_preserve_apd_nodes"] and row["all_minimizers_satisfy_boundary_targets"] for row in target_rows),
        "full_moment_uniqueness_is_fixed_support_conditional": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2582_payload = load_json(SOURCE_FILES["P2582_LOW_ORDER_MOMENT_MATCHING"])
    p2583_payload = load_json(SOURCE_FILES["P2583_FINITE_MOMENT_PREFIX_LADDER"])
    p2582 = theorem(p2582_payload, "apd_low_order_moment_matched_measure_nonuniqueness_audit")
    p2583 = theorem(p2583_payload, "apd_finite_moment_prefix_measure_ladder_audit")
    rows = apd_rows(p2416_payload)
    audit = full_moment_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2584_T1_apd_full_moment_finite_support_conditional_uniqueness_audit",
        "audited_chain": ["P2416/S1366", "P2582/S1532", "P2583/S1533"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "select the APD Gram measure by full moments on an already fixed finite support",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2582_low_order_moment_nonuniqueness_inherited": p2582.get("finite_low_order_moments_do_not_select_positive_measure") is True,
        "p2583_finite_moment_prefix_nonuniqueness_inherited": p2583.get("finite_moment_prefixes_do_not_select_positive_measure") is True,
        "apd_node_rows": rows,
        "apd_full_moment_finite_support_conditional_uniqueness_audit": audit,
        "full_moments_select_weights_only_after_support_is_fixed": audit["all_fixed_support_weights_conditionally_unique"],
        "finite_values_boundaries_and_full_support_moments_do_not_select_support": audit["all_targets_support_dependent"],
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Use full finite-support moments only as a conditional lemma: once the support is strict-sourced, Vandermonde moments recover the weights. P2584 shows that different strict-unsourced supports each have internally unique weights but still select different APD dynamics. The next honest step is to derive the APD measure support/density from strict nadsoliton dynamics, not merely solve weights on a chosen support."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2582_low_order_moment_nonuniqueness_inherited": theorem_export["p2582_low_order_moment_nonuniqueness_inherited"],
        "p2583_finite_moment_prefix_nonuniqueness_inherited": theorem_export["p2583_finite_moment_prefix_nonuniqueness_inherited"],
        "three_supports_audited": audit["support_count"] == 3,
        "three_targets_audited": audit["target_count"] == 3,
        "all_fixed_supports_have_vandermonde_rank_four": audit["all_fixed_supports_have_vandermonde_rank_four"],
        "all_fixed_support_weights_conditionally_unique": audit["all_fixed_support_weights_conditionally_unique"],
        "all_recovered_weights_positive": audit["all_recovered_weights_positive"],
        "all_targets_support_dependent": audit["all_targets_support_dependent"],
        "all_gram_metrics_positive_definite": audit["all_gram_metrics_positive_definite"],
        "all_minimizers_preserve_nodes_and_boundaries": audit["all_minimizers_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2584",
        "stage_id": "S1534",
        "status": "P2584_APD_FULL_MOMENT_FINITE_SUPPORT_CONDITIONAL_UNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_full_moment_finite_support_conditional_uniqueness_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2582_LOW_ORDER_MOMENT_MATCHING": sha256_json(p2582_payload),
                "P2583_FINITE_MOMENT_PREFIX_LADDER": sha256_json(p2583_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_full_moment_finite_support_conditional_uniqueness_audit"]["theorem_export"]
    audit = t["apd_full_moment_finite_support_conditional_uniqueness_audit"]
    lines = [
        "# P2584/S1534 APD full-moment finite-support conditional uniqueness audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Supports audited: `{audit['support_count']}`.",
        f"- Full moment orders: `{audit['full_moment_orders']}`.",
        f"- Fixed-support weights conditionally unique: `{audit['all_fixed_support_weights_conditionally_unique']}`.",
        f"- Targets support-dependent: `{audit['all_targets_support_dependent']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "Full moments on a fixed finite support recover the weights by a nonsingular Vandermonde system.  This is a useful conditional lemma, but not a strict source: changing the still-unsourced support changes the induced Gram metric and selected off-node APD dynamics while preserving the same APD node and boundary constraints.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No finite-support source, full-moment-law source, support-selection source, positive-measure source, L2-inner-product source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_full_moment_finite_support_conditional_uniqueness_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2584/S1534 APD full-moment finite-support conditional uniqueness guard

`P2584/S1534` records the positive conditional lemma after the finite-prefix obstructions: on a fixed four-point support, moments through order `3` recover the weights by a nonsingular Vandermonde system.  The same audit shows why this still is not a strict APD source: different unsourced supports each have conditionally unique positive weights but select different off-node APD dynamics while preserving finite APD nodes and endpoint slopes.
""".strip()
    lag_section = """
## P2584/S1534 APD full-moment finite-support conditional uniqueness Ltotal guard

`P2584/S1534` blocks a role-bearing APD Gram term in `L_total` from being justified by solving full finite-support moments after choosing the support.  The support or full measure density must itself be strict-sourced before the recovered weights can carry APD dynamics.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2584/S1534 APD full-moment finite-support conditional uniqueness guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2584/S1534 APD full-moment finite-support conditional uniqueness Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
