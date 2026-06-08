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
OUT = GEN / "p2583_s1533_apd_finite_moment_prefix_measure_ladder_audit.json"
MD = GEN / "p2583_s1533_apd_finite_moment_prefix_measure_ladder_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2581_APD_GRAM_MEASURE_DEPENDENCE": GEN / "p2581_s1531_apd_gram_measure_moment_dependence_audit.json",
    "P2582_APD_LOW_ORDER_MOMENT_MATCHING": GEN / "p2582_s1532_apd_low_order_moment_matched_measure_nonuniqueness_audit.json",
}
PREFIX_MAX_ORDERS = [0, 1, 2, 3]
LAMBDA_FRACTIONS = [-0.5, 0.0, 0.5]
NEGATIVE_EXPORT_FLAGS = [
    "apd_finite_moment_prefix_source_exported",
    "apd_truncated_moment_problem_source_exported",
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
        "new_packet": "P2583|S1533|APD finite moment prefix|finite moment prefix.*APD",
        "intended_research_nonduplication": "APD.*moment prefix ladder|moment prefix ladder.*APD|APD.*finite moment.*measure|finite moment.*APD.*measure|APD.*truncated moment|truncated moment.*APD|APD.*moment hierarchy|moment hierarchy.*APD|APD.*prefix moments|prefix moments.*APD",
        "apd_precursors": "P2416|S1366|P2581|S1531|P2582|S1532|APD.*Gram measure|APD.*moment-matched|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def moment_matrix(points: np.ndarray, max_order: int) -> np.ndarray:
    return np.vstack([points ** order for order in range(max_order + 1)])


def moment_null_vector(points: np.ndarray, max_order: int) -> np.ndarray:
    _, _, vh = np.linalg.svd(moment_matrix(points, max_order))
    vector = vh[-1]
    norm = np.linalg.norm(vector)
    if norm == 0.0:
        raise ValueError("moment matrix has no nonzero null vector")
    return vector / norm


def safe_lambda_radius(weights: np.ndarray, direction: np.ndarray) -> float:
    bounds = []
    for weight, delta in zip(weights, direction):
        if delta > 0.0:
            bounds.append(weight / delta)
        elif delta < 0.0:
            bounds.append(weight / -delta)
    if not bounds:
        raise ValueError("null direction cannot perturb weights")
    return 0.75 * min(bounds)


def prefix_measure_variants(max_order: int) -> dict[str, Any]:
    support_count = max(3, max_order + 2)
    points = np.linspace(0.25, 10.75, support_count)
    base_weights = np.ones(support_count)
    direction = moment_null_vector(points, max_order)
    radius = safe_lambda_radius(base_weights, direction)
    reference_moments = moment_matrix(points, max_order) @ base_weights
    variants = []
    for fraction in LAMBDA_FRACTIONS:
        lam = fraction * radius
        weights = base_weights + lam * direction
        moments = moment_matrix(points, max_order) @ weights
        variants.append({
            "name": f"prefix_{max_order}_lambda_{fraction:+.1f}".replace("+", "plus_").replace("-", "minus_").replace(".", "p"),
            "max_moment_order": max_order,
            "lambda_fraction": fraction,
            "lambda": float(lam),
            "points": [float(value) for value in points],
            "weights": [float(value) for value in weights],
            "matched_moments": [float(value) for value in moments],
            "max_abs_moment_deviation_from_reference": float(np.max(np.abs(moments - reference_moments))),
            "all_weights_positive": bool(np.min(weights) > 0.0),
        })
    return {
        "max_moment_order": max_order,
        "support_count": support_count,
        "support_points": [float(value) for value in points],
        "base_weights": [float(value) for value in base_weights],
        "moment_null_vector": [float(value) for value in direction],
        "safe_lambda_radius": float(radius),
        "reference_moments": [float(value) for value in reference_moments],
        "measure_variants": variants,
    }


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
        "max_moment_order": measure["max_moment_order"],
        "lambda_fraction": measure["lambda_fraction"],
        "weights": measure["weights"],
        "matched_moments": measure["matched_moments"],
        "max_abs_moment_deviation_from_reference": measure["max_abs_moment_deviation_from_reference"],
        "gram_metric_eigenvalues": [float(value) for value in eigenvalues],
        "gram_metric_is_positive_definite": bool(np.min(eigenvalues) > 0.0),
        "selected_coefficients": [float(value) for value in coefficients],
        "probe_values": [float(selected(point)) for point in PROBE_POINTS],
        "max_abs_node_residual": float(np.max(np.abs(node_residuals))),
        "max_abs_boundary_slope_residual": float(max(abs(left_residual), abs(right_residual))),
    }


def target_prefix_row(target: dict[str, Any], base: Any, basis: list[Any], matrix: np.ndarray, rows: list[dict[str, Any]], prefix: dict[str, Any]) -> dict[str, Any]:
    measure_rows = [measure_row(measure, base, basis, matrix, target, rows) for measure in prefix["measure_variants"]]
    middle_probe_values = {round(row["probe_values"][1], 12) for row in measure_rows}
    return {
        "target_name": target["name"],
        "max_moment_order": prefix["max_moment_order"],
        "measure_rows": measure_rows,
        "measure_count": len(measure_rows),
        "distinct_middle_probe_values_after_rounding_1e_minus_12": len(middle_probe_values),
        "all_measures_share_prefix_moments": all(row["max_abs_moment_deviation_from_reference"] <= 1.0e-10 for row in measure_rows),
        "all_measures_positive": all(min(row["weights"]) > 0.0 for row in measure_rows),
        "all_gram_metrics_positive_definite": all(row["gram_metric_is_positive_definite"] for row in measure_rows),
        "all_minimizers_preserve_apd_nodes": all(row["max_abs_node_residual"] <= 1.0e-6 for row in measure_rows),
        "all_minimizers_satisfy_boundary_targets": all(row["max_abs_boundary_slope_residual"] <= 1.0e-6 for row in measure_rows),
        "finite_prefix_moments_do_not_select_measure": len(middle_probe_values) > 1,
    }


def prefix_row(max_order: int, base: Any, basis: list[Any], matrix: np.ndarray, rows: list[dict[str, Any]]) -> dict[str, Any]:
    prefix = prefix_measure_variants(max_order)
    target_rows = [target_prefix_row(target, base, basis, matrix, rows, prefix) for target in BOUNDARY_TARGETS]
    return {
        **prefix,
        "target_rows": target_rows,
        "all_targets_prefix_nonunique": all(row["finite_prefix_moments_do_not_select_measure"] for row in target_rows),
        "all_measures_positive": all(measure["all_weights_positive"] for measure in prefix["measure_variants"]),
        "all_measures_share_prefix_moments": all(measure["max_abs_moment_deviation_from_reference"] <= 1.0e-10 for measure in prefix["measure_variants"]),
        "all_gram_metrics_positive_definite": all(row["all_gram_metrics_positive_definite"] for row in target_rows),
        "all_minimizers_preserve_nodes_and_boundaries": all(row["all_minimizers_preserve_apd_nodes"] and row["all_minimizers_satisfy_boundary_targets"] for row in target_rows),
    }


def moment_prefix_ladder_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    basis = reference_basis(vanish)
    matrix = boundary_matrix(basis)
    prefix_rows = [prefix_row(order, base, basis, matrix, rows) for order in PREFIX_MAX_ORDERS]
    return {
        "numpy_version": np.__version__,
        "prefix_max_orders": PREFIX_MAX_ORDERS,
        "lambda_fractions": LAMBDA_FRACTIONS,
        "probe_points": PROBE_POINTS,
        "prefix_rows": prefix_rows,
        "prefix_count": len(prefix_rows),
        "all_prefixes_have_positive_moment_matched_measures": all(row["all_measures_positive"] for row in prefix_rows),
        "all_prefixes_share_their_moment_prefix": all(row["all_measures_share_prefix_moments"] for row in prefix_rows),
        "all_prefixes_measure_nonunique": all(row["all_targets_prefix_nonunique"] for row in prefix_rows),
        "all_gram_metrics_positive_definite": all(row["all_gram_metrics_positive_definite"] for row in prefix_rows),
        "all_minimizers_preserve_nodes_and_boundaries": all(row["all_minimizers_preserve_nodes_and_boundaries"] for row in prefix_rows),
        "finite_moment_prefix_is_unsourced_measure_selector": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2581_payload = load_json(SOURCE_FILES["P2581_APD_GRAM_MEASURE_DEPENDENCE"])
    p2582_payload = load_json(SOURCE_FILES["P2582_APD_LOW_ORDER_MOMENT_MATCHING"])
    p2581 = theorem(p2581_payload, "apd_gram_measure_moment_dependence_audit")
    p2582 = theorem(p2582_payload, "apd_low_order_moment_matched_measure_nonuniqueness_audit")
    rows = apd_rows(p2416_payload)
    audit = moment_prefix_ladder_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2583_T1_apd_finite_moment_prefix_measure_ladder_audit",
        "audited_chain": ["P2416/S1366", "P2581/S1531", "P2582/S1532"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "select the APD Gram measure by a finite prefix of raw moments",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2581_measure_dependence_inherited": p2581.get("finite_apd_values_and_boundary_targets_do_not_select_positive_measure") is True,
        "p2582_low_order_moment_nonuniqueness_inherited": p2582.get("finite_low_order_moments_do_not_select_positive_measure") is True,
        "apd_node_rows": rows,
        "apd_finite_moment_prefix_measure_ladder_audit": audit,
        "finite_moment_prefixes_do_not_select_positive_measure": audit["all_prefixes_measure_nonunique"],
        "finite_moment_prefix_selector_is_conditional_not_strict_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not stop at any finite audited moment prefix for the APD measure. P2583 shows a ladder of positive measures sharing moment prefixes up through orders 0, 1, 2, and 3 while still selecting different APD off-node dynamics. The next honest step is to derive the full measure density, an infinite determinate moment law, or an equivalent strict APD kinetic source from nadsoliton dynamics."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2581_measure_dependence_inherited": theorem_export["p2581_measure_dependence_inherited"],
        "p2582_low_order_moment_nonuniqueness_inherited": theorem_export["p2582_low_order_moment_nonuniqueness_inherited"],
        "four_prefixes_audited": audit["prefix_count"] == 4,
        "all_prefixes_have_positive_moment_matched_measures": audit["all_prefixes_have_positive_moment_matched_measures"],
        "all_prefixes_share_their_moment_prefix": audit["all_prefixes_share_their_moment_prefix"],
        "all_prefixes_measure_nonunique": audit["all_prefixes_measure_nonunique"],
        "all_gram_metrics_positive_definite": audit["all_gram_metrics_positive_definite"],
        "all_minimizers_preserve_nodes_and_boundaries": audit["all_minimizers_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2583",
        "stage_id": "S1533",
        "status": "P2583_APD_FINITE_MOMENT_PREFIX_MEASURE_LADDER_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_finite_moment_prefix_measure_ladder_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2581_APD_GRAM_MEASURE_DEPENDENCE": sha256_json(p2581_payload),
                "P2582_APD_LOW_ORDER_MOMENT_MATCHING": sha256_json(p2582_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_finite_moment_prefix_measure_ladder_audit"]["theorem_export"]
    audit = t["apd_finite_moment_prefix_measure_ladder_audit"]
    lines = [
        "# P2583/S1533 APD finite moment-prefix measure ladder audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Moment prefixes audited: `{audit['prefix_max_orders']}`.",
        f"- All prefixes positive and moment-matched: `{audit['all_prefixes_have_positive_moment_matched_measures'] and audit['all_prefixes_share_their_moment_prefix']}`.",
        f"- All prefixes measure nonunique: `{audit['all_prefixes_measure_nonunique']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "Finite moment prefixes do not determine the APD Gram measure.  P2583 constructs a ladder of positive measures that share raw moment prefixes through orders 0, 1, 2, and 3; all constraints remain preserved, but off-node APD dynamics still varies across measures in every audited prefix.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No finite APD moment-prefix source, truncated-moment-problem source, positive-measure source, L2-inner-product source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_finite_moment_prefix_measure_ladder_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2583/S1533 APD finite moment-prefix measure ladder guard

`P2583/S1533` extends P2582 from one low-order prefix to a finite-prefix ladder.  For moment prefixes through orders `0`, `1`, `2`, and `3`, the audit constructs positive measures that share the audited raw moments and keep APD node/boundary constraints fixed, yet still select different off-node APD dynamics.  Thus no finite audited moment prefix is a strict APD measure source.
""".strip()
    lag_section = """
## P2583/S1533 APD finite moment-prefix measure ladder Ltotal guard

`P2583/S1533` blocks a role-bearing APD Gram term in `L_total` from being justified by any finite audited prefix of measure moments.  The strict action must provide a full density, determinate infinite moment law, or equivalent kinetic source before the APD inner product can carry dynamics.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2583/S1533 APD finite moment-prefix measure ladder guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2583/S1533 APD finite moment-prefix measure ladder Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
