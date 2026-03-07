#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "R8": load_json("fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs_summary.json"),
        "R9": load_json("fundamental_action_reconstruction/generated/r9_existing_kernel_feedback_host_to_control_pushforward_packet_for_kobs_summary.json"),
        "P10": load_json("fundamental_action_reconstruction/generated/p10_existing_kernel_feedback_to_kobs_rerun_after_explicit_current_pair_chain_summary.json"),
        "P13": load_json("fundamental_action_reconstruction/generated/p13_existing_kernel_feedback_factorization_rerun_after_host_to_control_pushforward_packet_summary.json"),
        "H8": load_json("fundamental_action_reconstruction/generated/h8_minimal_component_carrier_construction_spec_summary.json"),
        "H15": load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json"),
        "H33": load_json("fundamental_action_reconstruction/generated/h33_pair1_selector_target_justification_audit.json"),
        "H16": load_json("fundamental_action_reconstruction/generated/h16_operator_origin_partial_witness_audit_summary.json"),
    }

    checks_spec = [
        {
            "id": "r8_host_carrier_present",
            "actual": sources["R8"]["result"],
            "expected": "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "meaning": "the explicit operator-level host carrier is present",
        },
        {
            "id": "r9_typed_host_to_control_pushforward_present",
            "actual": sources["R9"]["result"],
            "expected": "typed_host_to_control_pushforward_present_but_control_scope_only_and_not_yet_reduced_to_the_explicit_current_pair_chain",
            "meaning": "the typed host-to-control pushforward is present",
        },
        {
            "id": "p10_full_block_still_present",
            "actual": sources["P10"]["status"],
            "expected": "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "meaning": "the explicit current-pair H3 block is still present on the route",
        },
        {
            "id": "p13_factorization_route_still_negative",
            "actual": sources["P13"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_TO_EXPLICIT_CHAIN_FACTORIZATION_ROUTE_AFTER_HOST_TO_CONTROL_PUSHFORWARD_PACKET",
            "meaning": "the direct factorization route remains noncomputable after the host-to-control pushforward packet",
        },
        {
            "id": "h8_explicit_chain_starts_on_pair1",
            "actual": sources["H8"]["pair"]["plane"],
            "expected": "V_1 = span{c1, s1}",
            "meaning": "the explicit H3 chain still starts on pair1 rather than on the full control carrier",
        },
        {
            "id": "h15_no_selector_sector_export_present",
            "actual": sources["H15"]["result"],
            "expected": "existing_feedback_not_identified_with_kobs",
            "meaning": "existing kernel feedback still has no selector-sector reduction/export of its own",
        },
        {
            "id": "h33_pair1_local_chart_only",
            "actual": sources["H33"]["result"],
            "expected": "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "meaning": "pair1 remains only a local chart and not a privileged selector target",
        },
        {
            "id": "h16_only_partial_witnesses_present",
            "actual": sources["H16"]["result"],
            "expected": "admissible_values_have_asymmetric_partial_witnesses_only",
            "meaning": "only partial operator-origin witnesses remain on the legacy side",
        },
    ]

    checks: list[dict[str, Any]] = []
    mismatches: list[str] = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    out = (
        ROOT
        / "generated"
        / "n16_current_existing_kernel_feedback_factorization_obstruction_after_host_to_control_pushforward_packet_theorem_summary.json"
    )

    if mismatches:
        summary = {
            "step": "N16",
            "status": "N16_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FACTORIZATION_ROUTE_FRONTIER_AFTER_HOST_TO_CONTROL_PUSHFORWARD_PACKET",
            "goal": "Check whether the current repo identifies existing kernel feedback with the explicit selector-facing H3 chain after host-to-control pushforward packetization.",
            "scope": "current_existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_R9_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected factorization-route frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_FACTORIZATION_ROUTE_FRONTIER_AFTER_R9_BEFORE_CLAIMING_N16",
        }
    else:
        summary = {
            "step": "N16",
            "status": "N16_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_FACTORIZATION_OBSTRUCTION_AFTER_HOST_TO_CONTROL_PUSHFORWARD_PACKET_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after materializing the typed host-to-control pushforward, the current repo still does not identify existing kernel feedback with the explicit selector-facing H3 chain.",
            "scope": "current_existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_R9_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "host_scope_operator_level_existing_kernel_feedback_carrier_present": True,
                "typed_host_to_control_pushforward_present": True,
                "full_current_pair_h3_block_present": True,
                "existing_kernel_feedback_identified_with_explicit_chain": False,
            },
            "missing_structure_classes": [
                "selector_sector_reduction_of_the_legacy_control_side_onto_pair1_or_an_equivalent_actual_target",
                "intertwiner_or_equality_witness_identifying_the_reduced_legacy_object_with_the_computed_current_pair_H3_block",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future factorization routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "IMPLEMENT_SELECTOR_SECTOR_REDUCTION_OF_THE_LEGACY_CONTROL_SIDE_AND_RERUN_P13",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
