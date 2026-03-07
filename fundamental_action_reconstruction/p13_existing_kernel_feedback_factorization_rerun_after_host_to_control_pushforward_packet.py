#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p13_existing_kernel_feedback_factorization_rerun_after_host_to_control_pushforward_packet.json"
OUT_SUMMARY = GENERATED / "p13_existing_kernel_feedback_factorization_rerun_after_host_to_control_pushforward_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p10 = load_json(
        "fundamental_action_reconstruction/generated/p10_existing_kernel_feedback_to_kobs_rerun_after_explicit_current_pair_chain_summary.json"
    )
    p12 = load_json(
        "fundamental_action_reconstruction/generated/p12_existing_kernel_feedback_factorization_rerun_after_host_carrier_packet.json"
    )
    r7 = load_json("fundamental_action_reconstruction/generated/r7_shared_frozen_kernel_provenance_packet_for_kobs_summary.json")
    r8 = load_json("fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs_summary.json")
    r9 = load_json("fundamental_action_reconstruction/generated/r9_existing_kernel_feedback_host_to_control_pushforward_packet_for_kobs_summary.json")
    h8 = load_json("fundamental_action_reconstruction/generated/h8_minimal_component_carrier_construction_spec_summary.json")
    h15 = load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json")
    h16 = load_json("fundamental_action_reconstruction/generated/h16_operator_origin_partial_witness_audit_summary.json")
    h33 = load_json("fundamental_action_reconstruction/generated/h33_pair1_selector_target_justification_audit.json")
    h34 = load_json("fundamental_action_reconstruction/generated/h34_basis_covariance_target_independence_audit.json")

    previous_missing = p12["remaining_missing_upstream_objects"]
    remaining_missing = [
        "selector_sector_reduction_of_the_legacy_control_side_onto_pair1_or_an_equivalent_actual_target",
        "intertwiner_or_equality_witness_identifying_the_reduced_legacy_object_with_the_computed_current_pair_H3_block",
    ]

    route_checks = [
        {
            "id": "p12_originally_missing_typed_projection",
            "pass": "typed_projection_or_pushforward_map_from_existing_kernel_feedback_into_the_explicit_H3_slot_chain_or_current_pair_block"
            in previous_missing,
            "expected": True,
            "actual": "typed_projection_or_pushforward_map_from_existing_kernel_feedback_into_the_explicit_H3_slot_chain_or_current_pair_block"
            in previous_missing,
            "meaning": "P12 indeed left the typed projection/pushforward as one remaining blocker",
        },
        {
            "id": "r7_shared_provenance_present",
            "pass": r7["result"] == "shared_frozen_kernel_provenance_present_but_not_operator_factorization",
            "expected": "shared_frozen_kernel_provenance_present_but_not_operator_factorization",
            "actual": r7["result"],
            "meaning": "shared frozen-kernel provenance remains present",
        },
        {
            "id": "r8_host_carrier_present",
            "pass": r8["result"]
            == "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "expected": "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "actual": r8["result"],
            "meaning": "the host-side legacy carrier remains present",
        },
        {
            "id": "r9_typed_host_to_control_pushforward_present",
            "pass": r9["result"]
            == "typed_host_to_control_pushforward_present_but_control_scope_only_and_not_yet_reduced_to_the_explicit_current_pair_chain",
            "expected": "typed_host_to_control_pushforward_present_but_control_scope_only_and_not_yet_reduced_to_the_explicit_current_pair_chain",
            "actual": r9["result"],
            "meaning": "the typed host-to-control pushforward is now present",
        },
        {
            "id": "h8_explicit_chain_starts_on_pair1",
            "pass": h8["pair"]["plane"] == "V_1 = span{c1, s1}",
            "expected": "V_1 = span{c1, s1}",
            "actual": h8["pair"]["plane"],
            "meaning": "the explicit H3 chain still starts on pair1 and not on the full control carrier",
        },
        {
            "id": "h33_pair1_local_chart_only",
            "pass": h33["result"]
            == "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "expected": "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "actual": h33["result"],
            "meaning": "pair1 is still only a local chart and not a privileged selector target",
        },
        {
            "id": "h34_no_covariance_argument",
            "pass": h34["status"] == "PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT",
            "expected": "PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT",
            "actual": h34["status"],
            "meaning": "no basis-covariant or target-independent reduction is present",
        },
        {
            "id": "h15_no_selector_sector_export_present",
            "pass": h15["result"] == "existing_feedback_not_identified_with_kobs",
            "expected": "existing_feedback_not_identified_with_kobs",
            "actual": h15["result"],
            "meaning": "existing kernel feedback still has no exported selector-sector reduction",
        },
        {
            "id": "h16_only_partial_witnesses_present",
            "pass": h16["result"] == "admissible_values_have_asymmetric_partial_witnesses_only",
            "expected": "admissible_values_have_asymmetric_partial_witnesses_only",
            "actual": h16["result"],
            "meaning": "the legacy operator-origin route still has only partial witnesses",
        },
        {
            "id": "p10_current_pair_block_still_present",
            "pass": p10["status"] == "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "expected": "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "actual": p10["status"],
            "meaning": "the explicit current-pair H3 block is still present on the route",
        },
    ]

    route_state = {
        "shared_frozen_kernel_provenance_present": True,
        "explicit_current_pair_h3_chain_present": True,
        "explicit_current_pair_h3_block_present": True,
        "explicit_operator_level_existing_kernel_feedback_carrier_present": True,
        "typed_host_to_control_pushforward_present": True,
        "legacy_control_carrier_present": True,
        "legacy_selector_sector_reduction_present": False,
        "intertwiner_or_equality_witness_present": False,
        "equivalence_factorization_map_present": False,
    }

    report = {
        "stage": "P13",
        "goal": "rerun_existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_host_to_control_pushforward_packet",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_TO_EXPLICIT_CHAIN_FACTORIZATION_ROUTE_AFTER_HOST_TO_CONTROL_PUSHFORWARD_PACKET",
        "reason": "the repo now has a typed host-to-control pushforward from the explicit legacy host carrier into the control-side mode carrier family, but the explicit H3 chain still starts on pair1 and the repo still lacks a selector-sector reduction from the legacy control side to pair1 or an equivalent actual target, plus an intertwiner/equality witness to the computed current-pair block",
        "lane": "existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_R9",
        "route_under_test": [
            "existing_kernel_feedback_inside_K_total",
            "shared_frozen_kernel_provenance_packet",
            "explicit_current_pair_H3_chain",
            "computed_current_pair_H3_block",
            "host_scope_operator_level_existing_kernel_feedback_carrier",
            "typed_host_to_control_pushforward",
            "equivalence_or_factorization_map",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "decomposition_of_P12_projection_blocker": {
            "from_P12": "typed_projection_or_pushforward_map_from_existing_kernel_feedback_into_the_explicit_H3_slot_chain_or_current_pair_block",
            "into_current_subobjects": [
                "typed_host_to_control_pushforward_from_the_legacy_host_carrier_into_M_control",
                "selector_sector_reduction_of_the_legacy_control_side_onto_pair1_or_an_equivalent_actual_target",
            ],
        },
        "resolved_from_P12": [
            "typed_host_to_control_pushforward_from_the_legacy_host_carrier_into_M_control"
        ],
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "P13_B1": "no selector-sector reduction currently exports the legacy control carrier onto pair1 or an equivalent actual target",
            "H15_B1": "existing kernel feedback still has no explicit selector-sector reduction or projected selector-block export in the current repository",
            "H16_B1": "both admissible operator-origin values still have only asymmetric partial witnesses and no equality/intertwiner witness to the computed current-pair block",
            "H33_B1": h33["frontier_text"],
            "H34_B1": h34["frontier_text"],
        },
        "required_next_step": "IMPLEMENT_SELECTOR_SECTOR_REDUCTION_OF_THE_LEGACY_CONTROL_SIDE_ONTO_PAIR1_OR_AN_EQUIVALENT_ACTUAL_TARGET_BEFORE_CLAIMING_ANY_EXISTING_KERNEL_FEEDBACK_FACTORIZATION",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P13",
        "status": report["status"],
        "reason": report["reason"],
        "resolved_from_P12": report["resolved_from_P12"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
