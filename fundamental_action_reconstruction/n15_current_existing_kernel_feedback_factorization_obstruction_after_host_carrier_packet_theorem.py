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
        "R7": load_json("fundamental_action_reconstruction/generated/r7_shared_frozen_kernel_provenance_packet_for_kobs_summary.json"),
        "R8": load_json("fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs_summary.json"),
        "P10": load_json("fundamental_action_reconstruction/generated/p10_existing_kernel_feedback_to_kobs_rerun_after_explicit_current_pair_chain_summary.json"),
        "P12": load_json("fundamental_action_reconstruction/generated/p12_existing_kernel_feedback_factorization_rerun_after_host_carrier_packet_summary.json"),
        "H15": load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json"),
        "H16": load_json("fundamental_action_reconstruction/generated/h16_operator_origin_partial_witness_audit_summary.json"),
    }

    checks_spec = [
        {
            "id": "r7_shared_provenance_present",
            "actual": sources["R7"]["result"],
            "expected": "shared_frozen_kernel_provenance_present_but_not_operator_factorization",
            "meaning": "shared frozen-kernel provenance is present",
        },
        {
            "id": "r8_host_carrier_present",
            "actual": sources["R8"]["result"],
            "expected": "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "meaning": "the explicit operator-level host carrier is now present",
        },
        {
            "id": "p10_full_block_still_present",
            "actual": sources["P10"]["status"],
            "expected": "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "meaning": "the explicit current-pair H3 block is still present on the route",
        },
        {
            "id": "p12_factorization_route_still_negative",
            "actual": sources["P12"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_TO_EXPLICIT_CHAIN_FACTORIZATION_ROUTE_AFTER_HOST_CARRIER_PACKET",
            "meaning": "the direct factorization route remains noncomputable after the host-carrier packet",
        },
        {
            "id": "h15_no_selector_sector_export_present",
            "actual": sources["H15"]["result"],
            "expected": "existing_feedback_not_identified_with_kobs",
            "meaning": "existing kernel feedback still has no selector-sector reduction/export of its own",
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
        / "n15_current_existing_kernel_feedback_factorization_obstruction_after_host_carrier_packet_theorem_summary.json"
    )

    if mismatches:
        summary = {
            "step": "N15",
            "status": "N15_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FACTORIZATION_ROUTE_FRONTIER_AFTER_HOST_CARRIER_PACKET",
            "goal": "Check whether the current repo identifies existing kernel feedback with the explicit selector-facing H3 chain after host-carrier packetization.",
            "scope": "current_existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_R8_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected factorization-route frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_FACTORIZATION_ROUTE_FRONTIER_AFTER_R8_BEFORE_CLAIMING_N15",
        }
    else:
        summary = {
            "step": "N15",
            "status": "N15_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_FACTORIZATION_OBSTRUCTION_AFTER_HOST_CARRIER_PACKET_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after materializing the host-scope operator-level existing-kernel-feedback carrier, the current repo still does not identify existing kernel feedback with the explicit selector-facing H3 chain.",
            "scope": "current_existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_R8_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "shared_frozen_kernel_provenance_present": True,
                "host_scope_operator_level_existing_kernel_feedback_carrier_present": True,
                "full_current_pair_h3_block_present": True,
                "existing_kernel_feedback_identified_with_explicit_chain": False,
            },
            "missing_structure_classes": [
                "typed_projection_or_pushforward_map_from_existing_kernel_feedback_into_the_explicit_H3_slot_chain_or_current_pair_block",
                "selector_sector_reduction_of_existing_kernel_feedback_onto_pair1_or_an_equivalent_actual_pair_target",
                "intertwiner_or_equality_witness_identifying_the_reduced_existing_feedback_object_with_the_computed_current_pair_H3_block",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future factorization routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "IMPLEMENT_TYPED_PROJECTION_FROM_THE_HOST_CARRIER_AND_RERUN_P12",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
