#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p14_existing_kernel_feedback_factorization_rerun_after_current_pair_chart_reduction_packet.json"
OUT_SUMMARY = GENERATED / "p14_existing_kernel_feedback_factorization_rerun_after_current_pair_chart_reduction_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p10 = load_json(
        "fundamental_action_reconstruction/generated/p10_existing_kernel_feedback_to_kobs_rerun_after_explicit_current_pair_chain_summary.json"
    )
    p13 = load_json(
        "fundamental_action_reconstruction/generated/p13_existing_kernel_feedback_factorization_rerun_after_host_to_control_pushforward_packet.json"
    )
    r8 = load_json("fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs_summary.json")
    r9 = load_json("fundamental_action_reconstruction/generated/r9_existing_kernel_feedback_host_to_control_pushforward_packet_for_kobs_summary.json")
    r10 = load_json("fundamental_action_reconstruction/generated/r10_legacy_control_to_current_pair_chart_reduction_packet_for_kobs_summary.json")
    h8 = load_json("fundamental_action_reconstruction/generated/h8_minimal_component_carrier_construction_spec_summary.json")
    h16 = load_json("fundamental_action_reconstruction/generated/h16_operator_origin_partial_witness_audit_summary.json")
    h33 = load_json("fundamental_action_reconstruction/generated/h33_pair1_selector_target_justification_audit.json")
    h34 = load_json("fundamental_action_reconstruction/generated/h34_basis_covariance_target_independence_audit.json")
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    previous_missing = p13["remaining_missing_upstream_objects"]
    remaining_missing = [
        "intertwiner_or_equality_witness_identifying_the_chart_reduced_legacy_object_with_the_computed_current_pair_H3_block"
    ]

    route_checks = [
        {
            "id": "p13_originally_missing_legacy_reduction",
            "pass": "selector_sector_reduction_of_the_legacy_control_side_onto_pair1_or_an_equivalent_actual_target"
            in previous_missing,
            "expected": True,
            "actual": "selector_sector_reduction_of_the_legacy_control_side_onto_pair1_or_an_equivalent_actual_target"
            in previous_missing,
            "meaning": "P13 indeed left the legacy reduction as one remaining blocker",
        },
        {
            "id": "r8_host_carrier_present",
            "pass": r8["result"]
            == "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "expected": "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "actual": r8["result"],
            "meaning": "the legacy host carrier remains present",
        },
        {
            "id": "r9_host_to_control_pushforward_present",
            "pass": r9["result"]
            == "typed_host_to_control_pushforward_present_but_control_scope_only_and_not_yet_reduced_to_the_explicit_current_pair_chain",
            "expected": "typed_host_to_control_pushforward_present_but_control_scope_only_and_not_yet_reduced_to_the_explicit_current_pair_chain",
            "actual": r9["result"],
            "meaning": "the legacy host-to-control pushforward remains present",
        },
        {
            "id": "r10_current_pair_chart_reduction_present",
            "pass": r10["result"]
            == "explicit_current_pair_chart_reduction_present_but_chart_scoped_only_and_not_a_strict_selector_target_justification",
            "expected": "explicit_current_pair_chart_reduction_present_but_chart_scoped_only_and_not_a_strict_selector_target_justification",
            "actual": r10["result"],
            "meaning": "the legacy control side now reaches the chosen explicit current-pair chart",
        },
        {
            "id": "h8_explicit_chain_starts_on_pair1",
            "pass": h8["pair"]["plane"] == "V_1 = span{c1, s1}",
            "expected": "V_1 = span{c1, s1}",
            "actual": h8["pair"]["plane"],
            "meaning": "the explicit H3 chain still starts on pair1",
        },
        {
            "id": "p10_current_pair_block_present",
            "pass": p10["status"] == "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "expected": "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "actual": p10["status"],
            "meaning": "the computed current-pair H3 block is still present on the route",
        },
        {
            "id": "c10_host_to_concrete_psi_block_identification_still_not_shown",
            "pass": c10["result"]["host_to_concrete_psi_block_identification"] == "not_shown",
            "expected": "not_shown",
            "actual": c10["result"]["host_to_concrete_psi_block_identification"],
            "meaning": "the repo still lacks coefficient-level identification from the legacy host to a concrete Psi block",
        },
        {
            "id": "h16_only_partial_witnesses_present",
            "pass": h16["result"] == "admissible_values_have_asymmetric_partial_witnesses_only",
            "expected": "admissible_values_have_asymmetric_partial_witnesses_only",
            "actual": h16["result"],
            "meaning": "the legacy operator-origin route still has only partial witnesses",
        },
        {
            "id": "h33_pair1_local_chart_only",
            "pass": h33["result"]
            == "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "expected": "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "actual": h33["result"],
            "meaning": "pair1 remains only a local chart",
        },
        {
            "id": "h34_no_covariance_argument",
            "pass": h34["status"] == "PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT",
            "expected": "PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT",
            "actual": h34["status"],
            "meaning": "no basis-covariant or target-independent promotion is present",
        },
    ]

    route_state = {
        "explicit_operator_level_existing_kernel_feedback_carrier_present": True,
        "typed_host_to_control_pushforward_present": True,
        "chosen_current_pair_chart_reduction_present": True,
        "computed_current_pair_h3_block_present": True,
        "chart_reduced_legacy_side_present": True,
        "intertwiner_or_equality_witness_present": False,
        "equivalence_factorization_map_present": False,
    }

    report = {
        "stage": "P14",
        "goal": "rerun_existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_current_pair_chart_reduction_packet",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_TO_EXPLICIT_CHAIN_FACTORIZATION_ROUTE_AFTER_CURRENT_PAIR_CHART_REDUCTION_PACKET",
        "reason": "the repo now carries the legacy route all the way to the chosen explicit current-pair chart used by the H3 chain, but it still lacks the operator-identification witness intertwining that chart-reduced legacy side with the computed current-pair H3 block",
        "lane": "existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_R10",
        "route_under_test": [
            "existing_kernel_feedback_inside_K_total",
            "shared_frozen_kernel_provenance_packet",
            "explicit_current_pair_H3_chain",
            "computed_current_pair_H3_block",
            "host_scope_operator_level_existing_kernel_feedback_carrier",
            "typed_host_to_control_pushforward",
            "current_pair_chart_reduction",
            "equivalence_or_factorization_map",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "resolved_from_P13": [
            "selector_sector_reduction_of_the_legacy_control_side_onto_the_chosen_explicit_current_pair_chart_pair1"
        ],
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "P14_B1": "no intertwiner or equality witness currently identifies the chart-reduced legacy object with the computed current-pair H3 block",
            "C10_B1": "no explicit coefficient-level or block-level identification between the certified legacy host operator and a concrete Psi-sector block",
            "H16_B1": "both admissible operator-origin values still have only asymmetric partial witnesses and no provenance-valid equality witness to the computed current-pair block",
            "H33_B1": h33["frontier_text"],
            "H34_B1": h34["frontier_text"],
        },
        "required_next_step": "IMPLEMENT_THE_SINGLE_REMAINING_INTERTWINER_OR_EQUALITY_WITNESS_BEFORE_CLAIMING_ANY_EXISTING_KERNEL_FEEDBACK_FACTORIZATION",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P14",
        "status": report["status"],
        "reason": report["reason"],
        "resolved_from_P13": report["resolved_from_P13"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
