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
from p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit import FULL_MOMENT_ORDERS, support_witness, target_row

GEN = ROOT / "generated"
OUT = GEN / "p2586_s1536_apd_mirror_symmetric_support_selector_nonuniqueness_audit.json"
MD = GEN / "p2586_s1536_apd_mirror_symmetric_support_selector_nonuniqueness_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2584_FULL_MOMENT_FIXED_SUPPORT": GEN / "p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit.json",
    "P2585_SUPPORT_GEOMETRY_NONUNIQUENESS": GEN / "p2585_s1535_apd_support_geometry_selector_nonuniqueness_audit.json",
}
MIRROR_CENTER = 5.5
MIRROR_PAIR_SUM = 2.0 * MIRROR_CENTER
MIRROR_SUPPORT_VARIANTS = [
    {"name": "mirror_internal_2p25_8p75", "points": [0.25, 2.25, 8.75, 10.75], "weights": [1.0, 1.0, 1.0, 1.0]},
    {"name": "mirror_internal_3p25_7p75", "points": [0.25, 3.25, 7.75, 10.75], "weights": [1.0, 1.0, 1.0, 1.0]},
    {"name": "mirror_internal_4p25_6p75", "points": [0.25, 4.25, 6.75, 10.75], "weights": [1.0, 1.0, 1.0, 1.0]},
]
NEGATIVE_EXPORT_FLAGS = [
    "apd_mirror_support_source_exported",
    "apd_reflection_symmetry_support_source_exported",
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
        "new_packet": "P2586|S1536|APD mirror support|mirror support.*APD",
        "intended_research_nonduplication": "APD.*symmetric support|symmetric support.*APD|APD.*reflection support|reflection support.*APD|APD.*support mirror|support mirror.*APD|APD.*pair-symmetric support|pair-symmetric support.*APD|APD.*support symmetry selector|support symmetry selector.*APD",
        "apd_precursors": "P2416|S1366|P2584|S1534|P2585|S1535|APD.*support geometry|APD.*finite support|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def mirror_support_witness(support: dict[str, Any]) -> dict[str, Any]:
    points = np.array(sorted(support["points"]), dtype=float)
    pair_sums = [float(points[0] + points[-1]), float(points[1] + points[-2])]
    return {
        **support_witness(support),
        "mirror_center": MIRROR_CENTER,
        "pair_sums": pair_sums,
        "support_centroid": float(np.mean(points)),
        "is_pair_symmetric_about_center": bool(all(abs(value - MIRROR_PAIR_SUM) <= 1.0e-12 for value in pair_sums)),
    }


def mirror_support_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    basis = reference_basis(vanish)
    matrix = boundary_matrix(basis)
    witnesses = [mirror_support_witness(support) for support in MIRROR_SUPPORT_VARIANTS]
    target_rows = [target_row(target, base, basis, matrix, rows, witnesses) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "mirror_center": MIRROR_CENTER,
        "mirror_pair_sum": MIRROR_PAIR_SUM,
        "full_moment_orders": FULL_MOMENT_ORDERS,
        "support_witnesses": witnesses,
        "target_rows": target_rows,
        "support_count": len(witnesses),
        "target_count": len(target_rows),
        "all_supports_pair_symmetric_about_center": all(witness["is_pair_symmetric_about_center"] for witness in witnesses),
        "all_fixed_support_weights_conditionally_unique": all(witness["full_moments_conditionally_select_weights_on_fixed_support"] for witness in witnesses),
        "all_fixed_supports_have_vandermonde_rank_four": all(witness["vandermonde_rank"] == 4 for witness in witnesses),
        "all_targets_mirror_support_nonunique": all(row["fixed_support_full_moments_do_not_select_support"] for row in target_rows),
        "all_gram_metrics_positive_definite": all(row["all_gram_metrics_positive_definite"] for row in target_rows),
        "all_minimizers_preserve_nodes_and_boundaries": all(row["all_minimizers_preserve_apd_nodes"] and row["all_minimizers_satisfy_boundary_targets"] for row in target_rows),
        "mirror_symmetry_is_unsourced_support_selector": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2584_payload = load_json(SOURCE_FILES["P2584_FULL_MOMENT_FIXED_SUPPORT"])
    p2585_payload = load_json(SOURCE_FILES["P2585_SUPPORT_GEOMETRY_NONUNIQUENESS"])
    p2584 = theorem(p2584_payload, "apd_full_moment_finite_support_conditional_uniqueness_audit")
    p2585 = theorem(p2585_payload, "apd_support_geometry_selector_nonuniqueness_audit")
    rows = apd_rows(p2416_payload)
    audit = mirror_support_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2586_T1_apd_mirror_symmetric_support_selector_nonuniqueness_audit",
        "audited_chain": ["P2416/S1366", "P2584/S1534", "P2585/S1535"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "select APD finite measure support by reflection/mirror pair symmetry plus fixed-support full moments",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2584_fixed_support_conditional_uniqueness_inherited": p2584.get("full_moments_select_weights_only_after_support_is_fixed") is True,
        "p2585_support_geometry_nonuniqueness_inherited": p2585.get("support_cardinality_endpoints_centroid_do_not_select_apd_support") is True,
        "apd_node_rows": rows,
        "apd_mirror_symmetric_support_selector_nonuniqueness_audit": audit,
        "mirror_pair_symmetry_does_not_select_apd_support": audit["all_targets_mirror_support_nonunique"],
        "fixed_support_full_moments_remain_conditional": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not promote mirror/reflection symmetry of finite support into an APD support source. P2586 keeps pair symmetry about the same center and still obtains distinct support choices with conditionally unique weights and support-dependent APD dynamics. The next honest step is to derive the actual support/density law from strict nadsoliton dynamics, not another geometric support preference."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2584_fixed_support_conditional_uniqueness_inherited": theorem_export["p2584_fixed_support_conditional_uniqueness_inherited"],
        "p2585_support_geometry_nonuniqueness_inherited": theorem_export["p2585_support_geometry_nonuniqueness_inherited"],
        "three_supports_audited": audit["support_count"] == 3,
        "three_targets_audited": audit["target_count"] == 3,
        "all_supports_pair_symmetric_about_center": audit["all_supports_pair_symmetric_about_center"],
        "all_fixed_support_weights_conditionally_unique": audit["all_fixed_support_weights_conditionally_unique"],
        "all_targets_mirror_support_nonunique": audit["all_targets_mirror_support_nonunique"],
        "all_gram_metrics_positive_definite": audit["all_gram_metrics_positive_definite"],
        "all_minimizers_preserve_nodes_and_boundaries": audit["all_minimizers_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2586",
        "stage_id": "S1536",
        "status": "P2586_APD_MIRROR_SYMMETRIC_SUPPORT_SELECTOR_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_mirror_symmetric_support_selector_nonuniqueness_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2584_FULL_MOMENT_FIXED_SUPPORT": sha256_json(p2584_payload),
                "P2585_SUPPORT_GEOMETRY_NONUNIQUENESS": sha256_json(p2585_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_mirror_symmetric_support_selector_nonuniqueness_audit"]["theorem_export"]
    audit = t["apd_mirror_symmetric_support_selector_nonuniqueness_audit"]
    lines = [
        "# P2586/S1536 APD mirror-symmetric support selector nonuniqueness audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Supports audited: `{audit['support_count']}`.",
        f"- Mirror center: `{audit['mirror_center']}`.",
        f"- All supports pair-symmetric: `{audit['all_supports_pair_symmetric_about_center']}`.",
        f"- Mirror support nonunique: `{audit['all_targets_mirror_support_nonunique']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "Mirror/pair symmetry is stronger than the centroid guard in P2585, but it still does not determine the APD support.  P2586 keeps reflection symmetry around the same center, recovers weights conditionally on each fixed support, and still obtains support-dependent off-node APD dynamics.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD mirror-support source, reflection-symmetry support source, support-selection source, finite-support source, positive-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_mirror_symmetric_support_selector_nonuniqueness_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2586/S1536 APD mirror-symmetric support selector nonuniqueness guard

`P2586/S1536` strengthens the P2585 support-geometry audit by imposing pairwise mirror symmetry about center `5.5`.  Three audited supports satisfy the same reflection law and each has conditionally unique positive weights from fixed-support full moments, but the selected off-node APD dynamics still changes with the support.  Thus mirror symmetry is not a strict APD support source.
""".strip()
    lag_section = """
## P2586/S1536 APD mirror-symmetric support selector nonuniqueness Ltotal guard

`P2586/S1536` blocks a role-bearing APD Gram term in `L_total` from being justified by finite-support reflection symmetry.  Symmetry of support is still a selector premise; strict nadsoliton dynamics must source the actual support/density law.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2586/S1536 APD mirror-symmetric support selector nonuniqueness guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2586/S1536 APD mirror-symmetric support selector nonuniqueness Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
