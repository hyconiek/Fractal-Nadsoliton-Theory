#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p12_existing_kernel_feedback_factorization_rerun_after_host_carrier_packet.json"
OUT_SUMMARY = GENERATED / "p12_existing_kernel_feedback_factorization_rerun_after_host_carrier_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p10 = load_json(
        "fundamental_action_reconstruction/generated/p10_existing_kernel_feedback_to_kobs_rerun_after_explicit_current_pair_chain_summary.json"
    )
    p11 = load_json(
        "fundamental_action_reconstruction/generated/p11_existing_kernel_feedback_to_explicit_h3_chain_factorization_map_probe.json"
    )
    r7 = load_json("fundamental_action_reconstruction/generated/r7_shared_frozen_kernel_provenance_packet_for_kobs_summary.json")
    r8 = load_json("fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs.json")
    r8_summary = load_json(
        "fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs_summary.json"
    )
    h14 = load_json("fundamental_action_reconstruction/generated/h14_existing_kernel_feedback_vs_new_kobs_separation_audit_summary.json")
    h15 = load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json")
    h16 = load_json("fundamental_action_reconstruction/generated/h16_operator_origin_partial_witness_audit_summary.json")

    previous_missing = p11["missing_upstream_objects"]
    remaining_missing = [
        "typed_projection_or_pushforward_map_from_existing_kernel_feedback_into_the_explicit_H3_slot_chain_or_current_pair_block",
        "selector_sector_reduction_of_existing_kernel_feedback_onto_pair1_or_an_equivalent_actual_pair_target",
        "intertwiner_or_equality_witness_identifying_the_reduced_existing_feedback_object_with_the_computed_current_pair_H3_block",
    ]

    route_checks = [
        {
            "id": "p11_originally_missing_explicit_legacy_carrier",
            "pass": "explicit_operator_level_existing_kernel_feedback_carrier_with_declared_basis_or_finite_state_space"
            in previous_missing,
            "expected": True,
            "actual": "explicit_operator_level_existing_kernel_feedback_carrier_with_declared_basis_or_finite_state_space"
            in previous_missing,
            "meaning": "P11 indeed localized the explicit operator-level legacy carrier as a missing subobject",
        },
        {
            "id": "r8_host_carrier_packet_present",
            "pass": r8_summary["result"]
            == "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "expected": "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "actual": r8_summary["result"],
            "meaning": "R8 now provides a host-scope operator-level existing-kernel-feedback carrier packet",
        },
        {
            "id": "r8_host_carrier_has_declared_finite_state_space",
            "pass": int(r8["declared_host_carrier"]["finite_state_space_size"]) == 12,
            "expected": 12,
            "actual": int(r8["declared_host_carrier"]["finite_state_space_size"]),
            "meaning": "the new host carrier packet really has a declared finite 12-state carrier",
        },
        {
            "id": "r7_shared_provenance_present",
            "pass": r7["result"] == "shared_frozen_kernel_provenance_present_but_not_operator_factorization",
            "expected": "shared_frozen_kernel_provenance_present_but_not_operator_factorization",
            "actual": r7["result"],
            "meaning": "shared frozen-kernel provenance remains present",
        },
        {
            "id": "p10_current_pair_block_still_present",
            "pass": p10["status"] == "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "expected": "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "actual": p10["status"],
            "meaning": "the explicit current-pair H3 block is still present on the route",
        },
        {
            "id": "h14_no_equivalence_map_present",
            "pass": h14["status"] == "PASS_PARTIAL_EXISTING_FEEDBACK_RECOGNIZED_BUT_NOT_IDENTIFIED_WITH_KOBS",
            "expected": "PASS_PARTIAL_EXISTING_FEEDBACK_RECOGNIZED_BUT_NOT_IDENTIFIED_WITH_KOBS",
            "actual": h14["status"],
            "meaning": "existing kernel feedback remains recognized but not identified with K_obs",
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
    ]

    route_state = {
        "shared_frozen_kernel_provenance_present": True,
        "existing_internal_feedback_parameter_packet_present": True,
        "explicit_current_pair_h3_chain_present": True,
        "explicit_current_pair_h3_block_present": True,
        "explicit_operator_level_existing_kernel_feedback_carrier_present": True,
        "explicit_operator_level_existing_kernel_feedback_carrier_scope": "host_scope_only",
        "typed_projection_from_existing_feedback_to_explicit_chain_present": False,
        "existing_feedback_selector_sector_reduction_present": False,
        "intertwiner_or_equality_witness_present": False,
        "equivalence_factorization_map_present": False,
    }

    report = {
        "stage": "P12",
        "goal": "rerun_existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_host_carrier_packet",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_TO_EXPLICIT_CHAIN_FACTORIZATION_ROUTE_AFTER_HOST_CARRIER_PACKET",
        "reason": "the repo now has a real host-scope operator-level existing-kernel-feedback carrier on a declared finite 12-state host carrier, but it still lacks a typed projection into the explicit H3 chain, a selector-sector reduction on the legacy side, and an intertwiner/equality witness identifying the reduced legacy object with the computed current-pair block",
        "lane": "existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_R8",
        "route_under_test": [
            "existing_kernel_feedback_inside_K_total",
            "shared_frozen_kernel_provenance_packet",
            "explicit_current_pair_H3_chain",
            "computed_current_pair_H3_block",
            "host_scope_operator_level_existing_kernel_feedback_carrier",
            "equivalence_or_factorization_map",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "resolved_from_P11": [
            "explicit_operator_level_existing_kernel_feedback_carrier_with_declared_basis_or_finite_state_space"
        ],
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R8_B1": "a host-scope operator-level existing-kernel-feedback carrier is now present, but it is not yet projected into the explicit H3 chain",
            "P12_B1": "no typed projection or pushforward currently exports the host-level existing-kernel-feedback carrier into the explicit H3 slot chain or current-pair block",
            "H15_B1": "existing kernel feedback still has no explicit selector-sector reduction or projected selector-block export in the current repository",
            "H16_B1": "both admissible operator-origin values still have only asymmetric partial witnesses and no equality/intertwiner witness to the computed current-pair block",
        },
        "required_next_step": "IMPLEMENT_TYPED_PROJECTION_OR_PUSHFORWARD_FROM_THE_HOST_CARRIER_INTO_THE_EXPLICIT_H3_CHAIN_BEFORE_CLAIMING_ANY_EXISTING_KERNEL_FEEDBACK_FACTORIZATION",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P12",
        "status": report["status"],
        "reason": report["reason"],
        "resolved_from_P11": report["resolved_from_P11"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
