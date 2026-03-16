#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p16_existing_kernel_feedback_legacy_chart_reduced_operator_export_probe.json"
OUT_SUMMARY = GENERATED / "p16_existing_kernel_feedback_legacy_chart_reduced_operator_export_probe_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p15 = load_json(
        "fundamental_action_reconstruction/generated/p15_existing_kernel_feedback_intertwiner_equality_witness_probe.json"
    )
    r8 = load_json(
        "fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs_summary.json"
    )
    r9 = load_json(
        "fundamental_action_reconstruction/generated/r9_existing_kernel_feedback_host_to_control_pushforward_packet_for_kobs_summary.json"
    )
    r10 = load_json(
        "fundamental_action_reconstruction/generated/r10_legacy_control_to_current_pair_chart_reduction_packet_for_kobs_summary.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")
    c11 = load_json("fundamental_action_reconstruction/generated/c11_psi_sector_block_extraction_audit_summary.json")
    c12 = load_json("fundamental_action_reconstruction/generated/c12_minimal_psi_block_extraction_packet_summary.json")
    c15 = load_json("fundamental_action_reconstruction/generated/c15_control_only_pullback_submatrix_packet_summary.json")
    c20 = load_json("fundamental_action_reconstruction/generated/c20_finite_materialization_recipe_audit_summary.json")

    p15_missing = p15["remaining_missing_upstream_objects"]
    remaining_missing = [
        "host_to_concrete_Psi_sector_quadratic_block_identification_for_the_existing_kernel_feedback_host_operator"
    ]

    route_checks = [
        {
            "id": "p15_originally_missing_legacy_chart_reduced_operator_object",
            "pass": "explicit_coefficient_filled_legacy_chart_reduced_operator_object_on_the_chosen_current_pair_chart_pair1_or_equivalent_actual_target"
            in p15_missing,
            "expected": True,
            "actual": "explicit_coefficient_filled_legacy_chart_reduced_operator_object_on_the_chosen_current_pair_chart_pair1_or_equivalent_actual_target"
            in p15_missing,
            "meaning": "P15 indeed localized the first remaining blocker as the missing coefficient-filled legacy chart-reduced operator object",
        },
        {
            "id": "r8_host_operator_carrier_present",
            "pass": r8["result"]
            == "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "expected": "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "actual": r8["result"],
            "meaning": "the legacy host operator carrier is present at host scope",
        },
        {
            "id": "r9_typed_host_to_control_pushforward_present",
            "pass": r9["result"]
            == "typed_host_to_control_pushforward_present_but_control_scope_only_and_not_yet_reduced_to_the_explicit_current_pair_chain",
            "expected": "typed_host_to_control_pushforward_present_but_control_scope_only_and_not_yet_reduced_to_the_explicit_current_pair_chain",
            "actual": r9["result"],
            "meaning": "the typed host-to-control pushforward is present",
        },
        {
            "id": "r10_current_pair_chart_reduction_present",
            "pass": r10["result"]
            == "explicit_current_pair_chart_reduction_present_but_chart_scoped_only_and_not_a_strict_selector_target_justification",
            "expected": "explicit_current_pair_chart_reduction_present_but_chart_scoped_only_and_not_a_strict_selector_target_justification",
            "actual": r10["result"],
            "meaning": "the chosen current-pair chart reduction is present",
        },
        {
            "id": "c10_host_to_concrete_psi_block_identification_not_shown",
            "pass": c10["result"]["host_to_concrete_psi_block_identification"] == "not_shown",
            "expected": "not_shown",
            "actual": c10["result"]["host_to_concrete_psi_block_identification"],
            "meaning": "the legacy host is still not identified with a concrete Psi-sector block",
        },
        {
            "id": "c15_coefficient_filled_h_psi_psi_present_via_r12",
            "pass": c15["result"]["coefficient_filled_H_PsiPsi_present"] == "yes",
            "expected": "yes",
            "actual": c15["result"]["coefficient_filled_H_PsiPsi_present"],
            "meaning": "a coefficient-filled canonical Psi x Psi block H_PsiPsi is exported (declared scope; see R12 pointer in C15)",
        },
        {
            "id": "c15_coefficient_filled_m_control_present_via_p476",
            "pass": c15["result"]["coefficient_filled_M_control_present"] == "yes",
            "expected": "yes",
            "actual": c15["result"]["coefficient_filled_M_control_present"],
            "meaning": "a coefficient-filled declared control pullback M_control is exported (declared scope; see P476 pointer in C15)",
        },
        {
            "id": "c15_orientation_slice_restriction_still_absent",
            "pass": c15["result"]["orientation_slice_restriction_present"] == "not_shown",
            "expected": "not_shown",
            "actual": c15["result"]["orientation_slice_restriction_present"],
            "meaning": "no declared restriction M_control -> orientation slice is exported at C15 level",
        },
    ]

    route_state = {
        "legacy_host_operator_carrier_present": True,
        "typed_host_to_control_pushforward_present": True,
        "chosen_current_pair_chart_reduction_present": True,
        "formal_control_pullback_formula_present": True,
        "host_to_concrete_psi_block_identification_present": False,
        "coefficient_filled_psi_sector_block_export_present": True,
        "coefficient_filled_control_pullback_present": True,
        "coefficient_filled_legacy_chart_reduced_operator_present": False,
    }

    report = {
        "stage": "P16",
        "goal": "compute_or_fail_existing_kernel_feedback_legacy_chart_reduced_operator_export_on_pair1_after_P15",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_LEGACY_CHART_REDUCED_OPERATOR_EXPORT_ROUTE",
        "reason": "the repo now contains the host carrier, the typed host-to-control pushforward, the formal control pullback formula, the chosen current-pair chart reduction, and (in declared scope) both a coefficient-filled canonical Psi x Psi block H_PsiPsi (R12) and a coefficient-filled declared control pullback M_control (P476), but it still exports no host-to-concrete Psi-block identification matching the existing-feedback host operator to that canonical carrier (C10_B1), so no strict existing-feedback legacy chart-reduced operator object on pair1 can be promoted from these declared artifacts",
        "lane": "existing_kernel_feedback_legacy_chart_reduced_operator_export_route_after_P15",
        "route_under_test": [
            "existing_kernel_feedback_host_operator",
            "host_to_concrete_Psi_sector_block_identification",
            "coefficient_filled_Psi_sector_block_export",
            "formal_control_pullback_M_control",
            "coefficient_filled_control_pullback_M_control",
            "chosen_current_pair_chart_reduction_to_pair1",
            "coefficient_filled_legacy_chart_reduced_operator_object_on_pair1",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R8_host_scope_existing_feedback_carrier",
            "R9_typed_host_to_control_pushforward",
            "C15_formal_control_pullback_formula",
            "R10_current_pair_chart_reduction_map",
            "R12_coefficient_filled_canonical_H_PsiPsi (declared support)",
            "P476_coefficient_filled_M_control (declared control pullback)",
        ],
        "decomposition_of_P15_missing_object": {
            "from_P15": "explicit_coefficient_filled_legacy_chart_reduced_operator_object_on_the_chosen_current_pair_chart_pair1_or_equivalent_actual_target",
            "into_current_blockers": remaining_missing,
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "C10_B1": c10["residual_blockers"]["C10_B1"],
            "C15_B2": c15["residual_blockers"]["C15_B2"],
        },
        "required_next_step": "EITHER_EXPORT_A_HOST_TO_CANONICAL_PSI_BLOCK_MATCHING_WITNESS_IDENTIFYING_QW2186_WITH_THE_CANONICAL_CARRIER_OR_KEEP_THE_LEGACY_OPERATOR_EXPORT_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P16",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P15_missing_object": report["decomposition_of_P15_missing_object"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
