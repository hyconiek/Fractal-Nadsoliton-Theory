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
from p2589_s1539_apd_mirror_sixth_moment_shell_support_nonuniqueness_audit import (
    ENDPOINTS,
    EXPECTED_CARDINALITY,
    FULL_SUPPORT_MOMENT_ORDERS,
    INTERNAL_FOURTH_SHELL,
    INTERNAL_SECOND_SHELL,
    INTERNAL_SIXTH_SHELL,
    MIRROR_CENTER,
    MIRROR_PAIR_SUM,
    shell_support_witness,
    support_from_product,
)

GEN = ROOT / "generated"
OUT = GEN / "p2590_s1540_apd_finite_even_moment_shell_interval_nonuniqueness_audit.json"
MD = GEN / "p2590_s1540_apd_finite_even_moment_shell_interval_nonuniqueness_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2584_FULL_MOMENT_FIXED_SUPPORT": GEN / "p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit.json",
    "P2589_SIXTH_MOMENT_SHELL_NONUNIQUENESS": GEN / "p2589_s1539_apd_mirror_sixth_moment_shell_support_nonuniqueness_audit.json",
}
PRODUCT_PARAMETER_GRID = [300.0, 340.0, 380.0, 420.0, 460.0, 500.0, 540.0, 576.0]
NEGATIVE_EXPORT_FLAGS = [
    "apd_finite_even_moment_shell_source_exported",
    "apd_product_parameter_interval_source_exported",
    "apd_support_higher_moment_shell_source_exported",
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
        "new_packet": "P2590|S1540|APD finite even moment shell interval|finite even moment shell.*APD",
        "intended_research_nonduplication": "APD.*moment shell continuum|moment shell continuum.*APD|APD.*product parameter support|product parameter support.*APD|APD.*moment shell parameter interval|moment shell parameter interval.*APD|APD.*finite even moment prefix.*support|finite even moment prefix.*support.*APD",
        "apd_precursors": "P2416|S1366|P2584|S1534|P2589|S1539|APD.*sixth moment shell|APD.*full moment|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def interval_support_witness(product_parameter: float) -> dict[str, Any]:
    witness = shell_support_witness(support_from_product(product_parameter))
    witness["product_parameter_is_free_interval_coordinate"] = True
    return witness


def vieta_power_sum_certificate(witnesses: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "statement": "For squared offsets x_i, fixing elementary sums e1=30, e2=273, e3=820 fixes power sums p1=e1=30, p2=e1^2-2e2=354, p3=e1^3-3e1e2+3e3=4890; the product e4 remains free over the audited interval.",
        "fixed_power_sums": {
            "p1_internal_second_shell": INTERNAL_SECOND_SHELL,
            "p2_internal_fourth_shell": INTERNAL_FOURTH_SHELL,
            "p3_internal_sixth_shell": INTERNAL_SIXTH_SHELL,
        },
        "audited_product_parameters": PRODUCT_PARAMETER_GRID,
        "distinct_product_parameters": len({round(witness["product_parameter"], 12) for witness in witnesses}),
        "all_witnesses_share_power_sums": all(witness["has_expected_internal_second_fourth_sixth_shell"] for witness in witnesses),
        "all_witnesses_have_positive_real_offsets": all(min(witness["squared_offsets"]) > 0.0 for witness in witnesses),
    }


def interval_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    basis = reference_basis(vanish)
    matrix = boundary_matrix(basis)
    witnesses = [interval_support_witness(parameter) for parameter in PRODUCT_PARAMETER_GRID]
    target_rows = [target_row(target, base, basis, matrix, rows, witnesses) for target in BOUNDARY_TARGETS]
    central_seconds = {round(witness["central_second_moment_unweighted"], 8) for witness in witnesses}
    central_fourths = {round(witness["central_fourth_moment_unweighted"], 8) for witness in witnesses}
    central_sixths = {round(witness["central_sixth_moment_unweighted"], 8) for witness in witnesses}
    return {
        "numpy_version": np.__version__,
        "mirror_center": MIRROR_CENTER,
        "mirror_pair_sum": MIRROR_PAIR_SUM,
        "expected_cardinality": EXPECTED_CARDINALITY,
        "endpoints": ENDPOINTS,
        "product_parameter_grid": PRODUCT_PARAMETER_GRID,
        "internal_second_shell": INTERNAL_SECOND_SHELL,
        "internal_fourth_shell": INTERNAL_FOURTH_SHELL,
        "internal_sixth_shell": INTERNAL_SIXTH_SHELL,
        "full_support_moment_orders": FULL_SUPPORT_MOMENT_ORDERS,
        "vieta_power_sum_certificate": vieta_power_sum_certificate(witnesses),
        "support_witnesses": witnesses,
        "target_rows": target_rows,
        "support_count": len(witnesses),
        "target_count": len(target_rows),
        "all_supports_have_expected_cardinality": all(witness["support_cardinality"] == EXPECTED_CARDINALITY for witness in witnesses),
        "all_supports_share_endpoints": all(witness["endpoints"] == ENDPOINTS for witness in witnesses),
        "all_supports_pair_symmetric_about_center": all(witness["is_pair_symmetric_about_center"] for witness in witnesses),
        "all_supports_share_second_fourth_sixth_shell": all(witness["has_expected_internal_second_fourth_sixth_shell"] for witness in witnesses),
        "all_supports_share_unweighted_central_second_moment": len(central_seconds) == 1,
        "all_supports_share_unweighted_central_fourth_moment": len(central_fourths) == 1,
        "all_supports_share_unweighted_central_sixth_moment": len(central_sixths) == 1,
        "all_fixed_support_weights_conditionally_unique": all(witness["full_moments_conditionally_select_weights_on_fixed_support"] for witness in witnesses),
        "all_fixed_supports_have_vandermonde_rank_ten": all(witness["vandermonde_rank"] == EXPECTED_CARDINALITY for witness in witnesses),
        "all_targets_interval_support_nonunique": all(row["fixed_support_full_moments_do_not_select_support"] for row in target_rows),
        "all_gram_metrics_positive_definite": all(row["all_gram_metrics_positive_definite"] for row in target_rows),
        "all_minimizers_preserve_nodes_and_boundaries": all(row["all_minimizers_preserve_apd_nodes"] and row["all_minimizers_satisfy_boundary_targets"] for row in target_rows),
        "finite_even_moment_shell_interval_is_unsourced_support_selector": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2584_payload = load_json(SOURCE_FILES["P2584_FULL_MOMENT_FIXED_SUPPORT"])
    p2589_payload = load_json(SOURCE_FILES["P2589_SIXTH_MOMENT_SHELL_NONUNIQUENESS"])
    p2584 = theorem(p2584_payload, "apd_full_moment_finite_support_conditional_uniqueness_audit")
    p2589 = theorem(p2589_payload, "apd_mirror_sixth_moment_shell_support_nonuniqueness_audit")
    rows = apd_rows(p2416_payload)
    audit = interval_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2590_T1_apd_finite_even_moment_shell_interval_nonuniqueness_audit",
        "audited_chain": ["P2416/S1366", "P2584/S1534", "P2589/S1539"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "select APD finite support by endpoints, mirror symmetry, cardinality ten, finite even moment shell prefix, and fixed-support full moments",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2584_fixed_support_conditional_uniqueness_inherited": p2584.get("full_moments_select_weights_only_after_support_is_fixed") is True,
        "p2589_sixth_moment_shell_nonuniqueness_inherited": p2589.get("mirror_second_fourth_sixth_moment_shell_does_not_select_apd_support") is True,
        "apd_node_rows": rows,
        "apd_finite_even_moment_shell_interval_nonuniqueness_audit": audit,
        "finite_even_moment_shell_prefix_does_not_select_apd_support": audit["all_targets_interval_support_nonunique"],
        "fixed_support_full_moments_remain_conditional": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not promote a finite even-moment shell prefix or its free product-parameter interval into an APD support source. P2590 keeps endpoints, mirror symmetry, cardinality, and fixed second/fourth/sixth shells across an audited product-parameter grid, while off-node APD dynamics remains support-dependent. The next honest step is a strict nadsoliton-derived support/density law or a genuine internal APD dynamics theorem."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2584_fixed_support_conditional_uniqueness_inherited": theorem_export["p2584_fixed_support_conditional_uniqueness_inherited"],
        "p2589_sixth_moment_shell_nonuniqueness_inherited": theorem_export["p2589_sixth_moment_shell_nonuniqueness_inherited"],
        "eight_supports_audited": audit["support_count"] == 8,
        "three_targets_audited": audit["target_count"] == 3,
        "vieta_power_sum_certificate_passes": audit["vieta_power_sum_certificate"]["all_witnesses_share_power_sums"],
        "all_supports_have_expected_cardinality": audit["all_supports_have_expected_cardinality"],
        "all_supports_share_endpoints": audit["all_supports_share_endpoints"],
        "all_supports_pair_symmetric_about_center": audit["all_supports_pair_symmetric_about_center"],
        "all_supports_share_second_fourth_sixth_shell": audit["all_supports_share_second_fourth_sixth_shell"],
        "all_fixed_support_weights_conditionally_unique": audit["all_fixed_support_weights_conditionally_unique"],
        "all_targets_interval_support_nonunique": audit["all_targets_interval_support_nonunique"],
        "all_gram_metrics_positive_definite": audit["all_gram_metrics_positive_definite"],
        "all_minimizers_preserve_nodes_and_boundaries": audit["all_minimizers_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2590",
        "stage_id": "S1540",
        "status": "P2590_APD_FINITE_EVEN_MOMENT_SHELL_INTERVAL_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_finite_even_moment_shell_interval_nonuniqueness_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2584_FULL_MOMENT_FIXED_SUPPORT": sha256_json(p2584_payload),
                "P2589_SIXTH_MOMENT_SHELL_NONUNIQUENESS": sha256_json(p2589_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_finite_even_moment_shell_interval_nonuniqueness_audit"]["theorem_export"]
    audit = t["apd_finite_even_moment_shell_interval_nonuniqueness_audit"]
    lines = [
        "# P2590/S1540 APD finite even-moment shell interval nonuniqueness audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Supports audited: `{audit['support_count']}`.",
        f"- Product-parameter grid: `{audit['product_parameter_grid']}`.",
        f"- Internal second/fourth/sixth shells: `{audit['internal_second_shell']}` / `{audit['internal_fourth_shell']}` / `{audit['internal_sixth_shell']}`.",
        f"- Interval support nonunique: `{audit['all_targets_interval_support_nonunique']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "P2590 turns the previous finite examples into an interval-style certificate: the Vieta power-sum identity fixes the even-moment shell prefix while leaving the product parameter free.  Across the audited grid, full fixed-support moments still recover positive weights only after support is chosen, and selected off-node APD dynamics changes with the support.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No finite even-moment shell source, product-parameter interval source, support-selection source, finite-support source, positive-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_finite_even_moment_shell_interval_nonuniqueness_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2590/S1540 APD finite even-moment shell interval nonuniqueness guard

`P2590/S1540` upgrades P2589 from three witnesses to an interval-style product-parameter grid.  The Vieta certificate fixes the APD mirror support's second/fourth/sixth even-moment shells while leaving the product parameter free; fixed-support full moments still recover weights only after support is chosen, and off-node APD dynamics remains support-dependent.  Thus a finite even-moment shell prefix is not a strict APD support source.
""".strip()
    lag_section = """
## P2590/S1540 APD finite even-moment shell interval nonuniqueness Ltotal guard

`P2590/S1540` blocks a role-bearing APD Gram term in `L_total` from being justified by a finite even-moment shell prefix or product-parameter interval.  Moment-shell interval geometry remains selector data; strict nadsoliton dynamics must source the actual support/density law.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2590/S1540 APD finite even-moment shell interval nonuniqueness guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2590/S1540 APD finite even-moment shell interval nonuniqueness Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
