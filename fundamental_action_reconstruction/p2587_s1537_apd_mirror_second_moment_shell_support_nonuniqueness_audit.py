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
from p2578_s1528_apd_augmented_boundary_basis_metric_dependence_audit import BOUNDARY_TARGETS, base_and_vanish
from p2580_s1530_apd_inner_product_basis_covariance_requirement_certificate import boundary_matrix
from p2581_s1531_apd_gram_measure_moment_dependence_audit import reference_basis
from p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit import target_row

GEN = ROOT / "generated"
OUT = GEN / "p2587_s1537_apd_mirror_second_moment_shell_support_nonuniqueness_audit.json"
MD = GEN / "p2587_s1537_apd_mirror_second_moment_shell_support_nonuniqueness_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2584_FULL_MOMENT_FIXED_SUPPORT": GEN / "p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit.json",
    "P2586_MIRROR_SUPPORT_NONUNIQUENESS": GEN / "p2586_s1536_apd_mirror_symmetric_support_selector_nonuniqueness_audit.json",
}
MIRROR_CENTER = 5.5
MIRROR_PAIR_SUM = 2.0 * MIRROR_CENTER
ENDPOINTS = [0.25, 10.75]
EXPECTED_CARDINALITY = 6
SHELL_RADIUS_SQUARED = 25.0
FULL_SUPPORT_MOMENT_ORDERS = [0, 1, 2, 3, 4, 5]
MIRROR_SECOND_MOMENT_SUPPORTS = [
    {
        "name": "six_point_shell_offsets_4_3",
        "offsets": [4.0, 3.0],
        "points": [0.25, 1.5, 2.5, 8.5, 9.5, 10.75],
        "weights": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    },
    {
        "name": "six_point_shell_offsets_24f5_7f5",
        "offsets": [24.0 / 5.0, 7.0 / 5.0],
        "points": [0.25, 0.7, 4.1, 6.9, 10.3, 10.75],
        "weights": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    },
    {
        "name": "six_point_shell_offsets_75f17_40f17",
        "offsets": [75.0 / 17.0, 40.0 / 17.0],
        "points": [0.25, 37.0 / 34.0, 107.0 / 34.0, 267.0 / 34.0, 337.0 / 34.0, 10.75],
        "weights": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    },
]
NEGATIVE_EXPORT_FLAGS = [
    "apd_mirror_second_moment_shell_source_exported",
    "apd_support_variance_shell_source_exported",
    "apd_support_selection_source_exported",
    "apd_finite_support_source_exported",
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
        "new_packet": "P2587|S1537|APD mirror second moment|mirror second moment.*APD",
        "intended_research_nonduplication": "APD.*second moment shell|second moment shell.*APD|APD.*support variance shell|support variance shell.*APD|APD.*central support moment|central support moment.*APD|APD.*mirror variance support|mirror variance support.*APD",
        "apd_precursors": "P2416|S1366|P2584|S1534|P2586|S1536|APD.*mirror.*support|APD.*full moment|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def vandermonde_moment_matrix(points: np.ndarray) -> np.ndarray:
    return np.vstack([points ** order for order in FULL_SUPPORT_MOMENT_ORDERS])


def full_support_moments(points: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return vandermonde_moment_matrix(points) @ weights


def shell_support_witness(support: dict[str, Any]) -> dict[str, Any]:
    points = np.array(sorted(support["points"]), dtype=float)
    weights = np.array(support["weights"], dtype=float)
    matrix = vandermonde_moment_matrix(points)
    moments = full_support_moments(points, weights)
    recovered = np.linalg.solve(matrix, moments)
    pair_sums = [float(points[index] + points[-1 - index]) for index in range(len(points) // 2)]
    declared_offsets = np.array(support["offsets"], dtype=float)
    internal_offset_square_sum = float(np.sum(declared_offsets ** 2))
    central_second_moment = float(np.sum((points - MIRROR_CENTER) ** 2))
    return {
        "support_name": support["name"],
        "points": [float(value) for value in points],
        "input_weights": [float(value) for value in weights],
        "declared_offsets": [float(value) for value in support["offsets"]],
        "support_cardinality": int(len(points)),
        "endpoints": [float(np.min(points)), float(np.max(points))],
        "mirror_center": MIRROR_CENTER,
        "pair_sums": pair_sums,
        "is_pair_symmetric_about_center": bool(all(abs(value - MIRROR_PAIR_SUM) <= 1.0e-12 for value in pair_sums)),
        "internal_offset_square_sum": internal_offset_square_sum,
        "has_expected_internal_second_moment_shell": bool(abs(internal_offset_square_sum - SHELL_RADIUS_SQUARED) <= 1.0e-10),
        "central_second_moment_unweighted": central_second_moment,
        "full_moment_orders": FULL_SUPPORT_MOMENT_ORDERS,
        "full_moments": [float(value) for value in moments],
        "vandermonde_determinant": float(np.linalg.det(matrix)),
        "vandermonde_rank": int(np.linalg.matrix_rank(matrix)),
        "recovered_weights_from_full_moments": [float(value) for value in recovered],
        "max_abs_recovered_weight_error": float(np.max(np.abs(recovered - weights))),
        "full_moments_conditionally_select_weights_on_fixed_support": bool(np.max(np.abs(recovered - weights)) <= 1.0e-8),
        "all_recovered_weights_positive": bool(np.min(recovered) > 0.0),
    }


def shell_support_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    basis = reference_basis(vanish)
    matrix = boundary_matrix(basis)
    witnesses = [shell_support_witness(support) for support in MIRROR_SECOND_MOMENT_SUPPORTS]
    target_rows = [target_row(target, base, basis, matrix, rows, witnesses) for target in BOUNDARY_TARGETS]
    central_seconds = {round(witness["central_second_moment_unweighted"], 10) for witness in witnesses}
    return {
        "numpy_version": np.__version__,
        "mirror_center": MIRROR_CENTER,
        "mirror_pair_sum": MIRROR_PAIR_SUM,
        "expected_cardinality": EXPECTED_CARDINALITY,
        "endpoints": ENDPOINTS,
        "internal_second_moment_shell_radius_squared": SHELL_RADIUS_SQUARED,
        "full_support_moment_orders": FULL_SUPPORT_MOMENT_ORDERS,
        "support_witnesses": witnesses,
        "target_rows": target_rows,
        "support_count": len(witnesses),
        "target_count": len(target_rows),
        "all_supports_have_expected_cardinality": all(witness["support_cardinality"] == EXPECTED_CARDINALITY for witness in witnesses),
        "all_supports_share_endpoints": all(witness["endpoints"] == ENDPOINTS for witness in witnesses),
        "all_supports_pair_symmetric_about_center": all(witness["is_pair_symmetric_about_center"] for witness in witnesses),
        "all_supports_share_second_moment_shell": all(witness["has_expected_internal_second_moment_shell"] for witness in witnesses),
        "all_supports_share_unweighted_central_second_moment": len(central_seconds) == 1,
        "all_fixed_support_weights_conditionally_unique": all(witness["full_moments_conditionally_select_weights_on_fixed_support"] for witness in witnesses),
        "all_fixed_supports_have_vandermonde_rank_six": all(witness["vandermonde_rank"] == 6 for witness in witnesses),
        "all_targets_shell_support_nonunique": all(row["fixed_support_full_moments_do_not_select_support"] for row in target_rows),
        "all_gram_metrics_positive_definite": all(row["all_gram_metrics_positive_definite"] for row in target_rows),
        "all_minimizers_preserve_nodes_and_boundaries": all(row["all_minimizers_preserve_apd_nodes"] and row["all_minimizers_satisfy_boundary_targets"] for row in target_rows),
        "second_moment_shell_is_unsourced_support_selector": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2584_payload = load_json(SOURCE_FILES["P2584_FULL_MOMENT_FIXED_SUPPORT"])
    p2586_payload = load_json(SOURCE_FILES["P2586_MIRROR_SUPPORT_NONUNIQUENESS"])
    p2584 = theorem(p2584_payload, "apd_full_moment_finite_support_conditional_uniqueness_audit")
    p2586 = theorem(p2586_payload, "apd_mirror_symmetric_support_selector_nonuniqueness_audit")
    rows = apd_rows(p2416_payload)
    audit = shell_support_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2587_T1_apd_mirror_second_moment_shell_support_nonuniqueness_audit",
        "audited_chain": ["P2416/S1366", "P2584/S1534", "P2586/S1536"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "select APD finite support by endpoints, mirror symmetry, cardinality six, fixed support second-moment shell, and fixed-support full moments",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2584_fixed_support_conditional_uniqueness_inherited": p2584.get("full_moments_select_weights_only_after_support_is_fixed") is True,
        "p2586_mirror_support_nonuniqueness_inherited": p2586.get("mirror_pair_symmetry_does_not_select_apd_support") is True,
        "apd_node_rows": rows,
        "apd_mirror_second_moment_shell_support_nonuniqueness_audit": audit,
        "mirror_second_moment_shell_does_not_select_apd_support": audit["all_targets_shell_support_nonunique"],
        "fixed_support_full_moments_remain_conditional": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not promote endpoints, mirror symmetry, cardinality, or a fixed second-moment support shell into an APD support source. P2587 keeps all of those constraints and still obtains distinct supports with conditionally unique weights and support-dependent APD dynamics. The next honest step is a strict nadsoliton-derived support/density law or a genuine internal APD dynamics theorem."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2584_fixed_support_conditional_uniqueness_inherited": theorem_export["p2584_fixed_support_conditional_uniqueness_inherited"],
        "p2586_mirror_support_nonuniqueness_inherited": theorem_export["p2586_mirror_support_nonuniqueness_inherited"],
        "three_supports_audited": audit["support_count"] == 3,
        "three_targets_audited": audit["target_count"] == 3,
        "all_supports_have_expected_cardinality": audit["all_supports_have_expected_cardinality"],
        "all_supports_share_endpoints": audit["all_supports_share_endpoints"],
        "all_supports_pair_symmetric_about_center": audit["all_supports_pair_symmetric_about_center"],
        "all_supports_share_second_moment_shell": audit["all_supports_share_second_moment_shell"],
        "all_fixed_support_weights_conditionally_unique": audit["all_fixed_support_weights_conditionally_unique"],
        "all_targets_shell_support_nonunique": audit["all_targets_shell_support_nonunique"],
        "all_gram_metrics_positive_definite": audit["all_gram_metrics_positive_definite"],
        "all_minimizers_preserve_nodes_and_boundaries": audit["all_minimizers_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2587",
        "stage_id": "S1537",
        "status": "P2587_APD_MIRROR_SECOND_MOMENT_SHELL_SUPPORT_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_mirror_second_moment_shell_support_nonuniqueness_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2584_FULL_MOMENT_FIXED_SUPPORT": sha256_json(p2584_payload),
                "P2586_MIRROR_SUPPORT_NONUNIQUENESS": sha256_json(p2586_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_mirror_second_moment_shell_support_nonuniqueness_audit"]["theorem_export"]
    audit = t["apd_mirror_second_moment_shell_support_nonuniqueness_audit"]
    lines = [
        "# P2587/S1537 APD mirror second-moment shell support nonuniqueness audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Supports audited: `{audit['support_count']}`.",
        f"- Cardinality: `{audit['expected_cardinality']}`.",
        f"- Mirror center: `{audit['mirror_center']}`.",
        f"- Internal second-moment shell: `{audit['internal_second_moment_shell_radius_squared']}`.",
        f"- Shell support nonunique: `{audit['all_targets_shell_support_nonunique']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "P2587 strengthens P2586 by adding a fixed second-moment shell on six-point mirror supports.  The supports share endpoints, cardinality, reflection symmetry, and unweighted central second moment; fixed-support full moments recover positive weights, but selected off-node APD dynamics still changes with the support.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No second-moment shell support source, variance-shell source, support-selection source, finite-support source, positive-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_mirror_second_moment_shell_support_nonuniqueness_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2587/S1537 APD mirror second-moment shell support nonuniqueness guard

`P2587/S1537` adds a stronger support-geometry constraint than P2586: six-point supports share endpoints, mirror symmetry about center `5.5`, cardinality `6`, and a common internal second-moment shell.  Fixed-support full moments still recover positive weights only after support is chosen, and the off-node APD dynamics remains support-dependent.  Therefore a mirror variance shell is not a strict APD support source.
""".strip()
    lag_section = """
## P2587/S1537 APD mirror second-moment shell support nonuniqueness Ltotal guard

`P2587/S1537` blocks a role-bearing APD Gram term in `L_total` from being justified by a finite support variance shell.  Endpoints, reflection symmetry, cardinality, and second support moment are still selector data; strict nadsoliton dynamics must source the actual support/density law.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2587/S1537 APD mirror second-moment shell support nonuniqueness guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2587/S1537 APD mirror second-moment shell support nonuniqueness Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
