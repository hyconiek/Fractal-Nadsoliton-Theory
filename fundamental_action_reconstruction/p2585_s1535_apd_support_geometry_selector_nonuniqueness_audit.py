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
from p2581_s1531_apd_gram_measure_moment_dependence_audit import reference_basis
from p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit import FULL_MOMENT_ORDERS, support_witness, target_row
from p2580_s1530_apd_inner_product_basis_covariance_requirement_certificate import boundary_matrix

GEN = ROOT / "generated"
OUT = GEN / "p2585_s1535_apd_support_geometry_selector_nonuniqueness_audit.json"
MD = GEN / "p2585_s1535_apd_support_geometry_selector_nonuniqueness_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2583_FINITE_MOMENT_PREFIX_LADDER": GEN / "p2583_s1533_apd_finite_moment_prefix_measure_ladder_audit.json",
    "P2584_FULL_MOMENT_FIXED_SUPPORT": GEN / "p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit.json",
}
SUPPORT_GEOMETRY_VARIANTS = [
    {"name": "balanced_internal_3p25_7p25", "points": [0.25, 3.25, 7.25, 10.75], "weights": [1.0, 1.0, 1.0, 1.0]},
    {"name": "balanced_internal_4p25_6p25", "points": [0.25, 4.25, 6.25, 10.75], "weights": [1.0, 1.0, 1.0, 1.0]},
    {"name": "balanced_internal_2p25_8p25", "points": [0.25, 2.25, 8.25, 10.75], "weights": [1.0, 1.0, 1.0, 1.0]},
]
EXPECTED_SUPPORT_CARDINALITY = 4
EXPECTED_LEFT_ENDPOINT = 0.25
EXPECTED_RIGHT_ENDPOINT = 10.75
EXPECTED_SUPPORT_CENTROID = 5.375
NEGATIVE_EXPORT_FLAGS = [
    "apd_support_geometry_source_exported",
    "apd_support_endpoint_centroid_source_exported",
    "apd_support_cardinality_source_exported",
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
        "new_packet": "P2585|S1535|APD support geometry|support geometry.*APD",
        "intended_research_nonduplication": "APD.*support centroid|support centroid.*APD|APD.*support endpoint|support endpoint.*APD|APD.*support cardinality|support cardinality.*APD|APD.*finite support geometry|finite support geometry.*APD|APD.*support selector.*centroid|support selector.*centroid.*APD",
        "apd_precursors": "P2416|S1366|P2583|S1533|P2584|S1534|APD.*finite support|APD.*full moment|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def support_geometry_witness(support: dict[str, Any]) -> dict[str, Any]:
    points = np.array(support["points"], dtype=float)
    return {
        **support_witness(support),
        "support_cardinality": int(len(points)),
        "left_endpoint": float(np.min(points)),
        "right_endpoint": float(np.max(points)),
        "support_centroid": float(np.mean(points)),
        "has_expected_cardinality_endpoint_centroid": bool(
            len(points) == EXPECTED_SUPPORT_CARDINALITY
            and np.min(points) == EXPECTED_LEFT_ENDPOINT
            and np.max(points) == EXPECTED_RIGHT_ENDPOINT
            and abs(np.mean(points) - EXPECTED_SUPPORT_CENTROID) <= 1.0e-12
        ),
    }


def support_geometry_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    basis = reference_basis(vanish)
    matrix = boundary_matrix(basis)
    witnesses = [support_geometry_witness(support) for support in SUPPORT_GEOMETRY_VARIANTS]
    target_rows = [target_row(target, base, basis, matrix, rows, witnesses) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "expected_support_cardinality": EXPECTED_SUPPORT_CARDINALITY,
        "expected_left_endpoint": EXPECTED_LEFT_ENDPOINT,
        "expected_right_endpoint": EXPECTED_RIGHT_ENDPOINT,
        "expected_support_centroid": EXPECTED_SUPPORT_CENTROID,
        "full_moment_orders": FULL_MOMENT_ORDERS,
        "support_witnesses": witnesses,
        "target_rows": target_rows,
        "support_count": len(witnesses),
        "target_count": len(target_rows),
        "all_supports_share_geometry_constraints": all(witness["has_expected_cardinality_endpoint_centroid"] for witness in witnesses),
        "all_fixed_support_weights_conditionally_unique": all(witness["full_moments_conditionally_select_weights_on_fixed_support"] for witness in witnesses),
        "all_fixed_supports_have_vandermonde_rank_four": all(witness["vandermonde_rank"] == 4 for witness in witnesses),
        "all_targets_support_geometry_nonunique": all(row["fixed_support_full_moments_do_not_select_support"] for row in target_rows),
        "all_gram_metrics_positive_definite": all(row["all_gram_metrics_positive_definite"] for row in target_rows),
        "all_minimizers_preserve_nodes_and_boundaries": all(row["all_minimizers_preserve_apd_nodes"] and row["all_minimizers_satisfy_boundary_targets"] for row in target_rows),
        "support_geometry_constraints_are_unsourced_selector": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2583_payload = load_json(SOURCE_FILES["P2583_FINITE_MOMENT_PREFIX_LADDER"])
    p2584_payload = load_json(SOURCE_FILES["P2584_FULL_MOMENT_FIXED_SUPPORT"])
    p2583 = theorem(p2583_payload, "apd_finite_moment_prefix_measure_ladder_audit")
    p2584 = theorem(p2584_payload, "apd_full_moment_finite_support_conditional_uniqueness_audit")
    rows = apd_rows(p2416_payload)
    audit = support_geometry_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2585_T1_apd_support_geometry_selector_nonuniqueness_audit",
        "audited_chain": ["P2416/S1366", "P2583/S1533", "P2584/S1534"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "select APD finite measure support by cardinality, endpoints, centroid, and fixed-support full moments",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2583_finite_moment_prefix_nonuniqueness_inherited": p2583.get("finite_moment_prefixes_do_not_select_positive_measure") is True,
        "p2584_fixed_support_conditional_uniqueness_inherited": p2584.get("full_moments_select_weights_only_after_support_is_fixed") is True,
        "apd_node_rows": rows,
        "apd_support_geometry_selector_nonuniqueness_audit": audit,
        "support_cardinality_endpoints_centroid_do_not_select_apd_support": audit["all_targets_support_geometry_nonunique"],
        "fixed_support_full_moments_remain_conditional": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not promote finite-support geometry constraints into an APD support source. P2585 fixes cardinality, endpoints, and centroid, and still obtains distinct support choices with conditionally unique Vandermonde weights that select different APD off-node dynamics. The next honest step is to derive the actual APD support/density law from strict nadsoliton dynamics."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2583_finite_moment_prefix_nonuniqueness_inherited": theorem_export["p2583_finite_moment_prefix_nonuniqueness_inherited"],
        "p2584_fixed_support_conditional_uniqueness_inherited": theorem_export["p2584_fixed_support_conditional_uniqueness_inherited"],
        "three_supports_audited": audit["support_count"] == 3,
        "three_targets_audited": audit["target_count"] == 3,
        "all_supports_share_geometry_constraints": audit["all_supports_share_geometry_constraints"],
        "all_fixed_support_weights_conditionally_unique": audit["all_fixed_support_weights_conditionally_unique"],
        "all_targets_support_geometry_nonunique": audit["all_targets_support_geometry_nonunique"],
        "all_gram_metrics_positive_definite": audit["all_gram_metrics_positive_definite"],
        "all_minimizers_preserve_nodes_and_boundaries": audit["all_minimizers_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2585",
        "stage_id": "S1535",
        "status": "P2585_APD_SUPPORT_GEOMETRY_SELECTOR_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_support_geometry_selector_nonuniqueness_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2583_FINITE_MOMENT_PREFIX_LADDER": sha256_json(p2583_payload),
                "P2584_FULL_MOMENT_FIXED_SUPPORT": sha256_json(p2584_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_support_geometry_selector_nonuniqueness_audit"]["theorem_export"]
    audit = t["apd_support_geometry_selector_nonuniqueness_audit"]
    lines = [
        "# P2585/S1535 APD support-geometry selector nonuniqueness audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Supports audited: `{audit['support_count']}`.",
        f"- Shared support constraints: cardinality `{audit['expected_support_cardinality']}`, endpoints `({audit['expected_left_endpoint']}, {audit['expected_right_endpoint']})`, centroid `{audit['expected_support_centroid']}`.",
        f"- Fixed-support weights conditionally unique: `{audit['all_fixed_support_weights_conditionally_unique']}`.",
        f"- Support geometry nonunique: `{audit['all_targets_support_geometry_nonunique']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "Cardinality, endpoints, centroid, and fixed-support full moments still do not determine the APD support.  P2585 keeps those support-geometry constraints fixed, recovers weights conditionally on each support, and still obtains support-dependent off-node APD dynamics.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD support-geometry source, endpoint/centroid source, support-cardinality source, finite-support source, positive-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_support_geometry_selector_nonuniqueness_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2585/S1535 APD support-geometry selector nonuniqueness guard

`P2585/S1535` tests whether P2584's remaining support choice can be sourced by simple finite-support geometry.  Three audited supports share cardinality `4`, endpoints `0.25` and `10.75`, and centroid `5.375`; on each support, full moments conditionally recover positive weights.  The selected off-node APD dynamics still changes with the support, so support geometry is not `strict_dynamical_source_for_A_P_D`.
""".strip()
    lag_section = """
## P2585/S1535 APD support-geometry selector nonuniqueness Ltotal guard

`P2585/S1535` blocks a role-bearing APD Gram term in `L_total` from being justified by support cardinality, endpoints, and centroid.  These geometry constraints do not source the APD support/density; strict nadsoliton dynamics must provide the actual support law before recovered weights can carry dynamics.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2585/S1535 APD support-geometry selector nonuniqueness guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2585/S1535 APD support-geometry selector nonuniqueness Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
