#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p24_existing_kernel_feedback_host_matching_witness_rerun_after_host_side_residual_absence_packet.json"
OUT_SUMMARY = GENERATED / "p24_existing_kernel_feedback_host_matching_witness_rerun_after_host_side_residual_absence_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p23 = load_json(
        "fundamental_action_reconstruction/generated/p23_existing_kernel_feedback_host_matching_witness_rerun_after_residual_diagonal_pullback_packet.json"
    )
    r17 = load_json(
        "fundamental_action_reconstruction/generated/r17_explicit_host_side_residual_diagonal_correction_absence_packet_for_host_matching_route.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    prior_gap = "explicit_zero_or_host_side_cancellation_witness_for_the_declared_control_pullback_of_the_residual_local_diagonal_sector"
    remaining_missing = [
        "explicit_zero_witness_for_the_declared_control_pullback_of_the_residual_local_diagonal_sector",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    route_checks = [
        {
            "id": "p23_zero_or_host_side_gap_was_still_missing",
            "actual": prior_gap in p23["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P23 still had the zero-or-host-side cancellation gap as a missing object",
        },
        {
            "id": "r17_host_side_absence_packet_present",
            "actual": r17["result"]["explicit_host_side_residual_diagonal_correction_absence_packet_present"],
            "expected": True,
            "meaning": "R17 adds the explicit host-side residual diagonal correction absence packet",
        },
        {
            "id": "r17_host_side_residual_correction_absent",
            "actual": r17["result"]["host_side_residual_diagonal_correction_present"],
            "expected": False,
            "meaning": "R17 shows that no host-side residual diagonal correction exists on the current route",
        },
        {
            "id": "r17_zero_witness_for_canonical_residual_declared_pullback_still_absent",
            "actual": r17["result"]["zero_witness_for_canonical_residual_declared_pullback_present"],
            "expected": False,
            "meaning": "R17 does not claim that the canonical residual declared pullback is zero",
        },
        {
            "id": "q2191_full_physical_uniqueness_still_open",
            "actual": q2191["flags"]["full_physical_uniqueness_closed"],
            "expected": False,
            "meaning": "QW-2191 still blocks full physical uniqueness",
        },
        {
            "id": "c10_host_identification_not_shown",
            "actual": c10["result"]["host_to_concrete_psi_block_identification"],
            "expected": "not_shown",
            "meaning": "host-to-concrete Psi-block identification remains absent",
        },
    ]

    for item in route_checks:
        item["pass"] = item["actual"] == item["expected"]

    route_state = {
        "partial_host_to_canonical_block_overlap_present": True,
        "kernel_coefficient_specialization_witness_present": True,
        "host_scalar_floor_embedding_packet_present": True,
        "declared_control_pullback_of_residual_local_diagonal_sector_present": True,
        "host_side_residual_diagonal_correction_present": False,
        "zero_witness_for_canonical_residual_declared_pullback_present": False,
        "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
    }

    report = {
        "stage": "P24",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_matching_witness_after_R17",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_WITNESS_ROUTE_AFTER_R17_HOST_SIDE_RESIDUAL_ABSENCE_PACKET",
        "reason": "the repo now shows that no host-side residual diagonal correction exists beyond K_total plus m0^2 I, so the only remaining diagonal-side blocker is a zero witness for the canonical residual declared pullback; the QW-2191 canonicalization boundary also remains open",
        "lane": "existing_kernel_feedback_host_matching_witness_route_after_R17",
        "route_under_test": [
            "partial_host_to_canonical_block_overlap_packet",
            "kernel_coefficient_specialization_witness",
            "host_scalar_floor_embedding_packet",
            "declared_control_pullback_of_residual_local_diagonal_sector",
            "host_side_residual_diagonal_correction_branch",
            "zero_witness_for_the_canonical_residual_declared_pullback",
            "full_physical_uniqueness_or_selector_relevant_canonicalization",
            "host_to_concrete_Psi_block_identification",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R17_explicit_host_side_residual_diagonal_correction_absence_packet",
            "R16_explicit_residual_local_diagonal_declared_control_pullback_packet",
            "R15_explicit_host_scalar_floor_embedding_packet",
            "R14_explicit_frozen_kernel_specialization_packet",
        ],
        "decomposition_of_P23_zero_or_host_side_gap": {
            "from_P23": prior_gap,
            "discharged_by_R17_as": "explicit_host_side_residual_diagonal_correction_absence_packet",
            "remaining_matching_gap": "explicit_zero_witness_for_the_declared_control_pullback_of_the_residual_local_diagonal_sector",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R17_B1": r17["frontier_text"],
            "QW2191_required_next_step": q2191["required_next_step"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_SELECTOR_RELEVANT_CANONICALIZATION_AND_ADD_AN_EXPLICIT_ZERO_WITNESS_FOR_THE_CANONICAL_RESIDUAL_DECLARED_PULLBACK_OR_KEEP_THE_HOST_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P24",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P23_zero_or_host_side_gap": report["decomposition_of_P23_zero_or_host_side_gap"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
