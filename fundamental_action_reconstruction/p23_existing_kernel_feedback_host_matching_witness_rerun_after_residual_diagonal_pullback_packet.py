#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p23_existing_kernel_feedback_host_matching_witness_rerun_after_residual_diagonal_pullback_packet.json"
OUT_SUMMARY = GENERATED / "p23_existing_kernel_feedback_host_matching_witness_rerun_after_residual_diagonal_pullback_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p22 = load_json(
        "fundamental_action_reconstruction/generated/p22_existing_kernel_feedback_host_matching_witness_rerun_after_diagonal_floor_embedding_packet.json"
    )
    r16 = load_json(
        "fundamental_action_reconstruction/generated/r16_explicit_residual_local_diagonal_declared_control_pullback_packet_for_host_matching_route.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    prior_gap = (
        "explicit_residual_local_diagonal_sector_equality_or_cancellation_witness_reducing_the_canonical_diagonal_sector_to_the_host_floor_m0_squared_I_or_to_a_declared_control_pullback_of_it"
    )
    remaining_missing = [
        "explicit_zero_or_host_side_cancellation_witness_for_the_declared_control_pullback_of_the_residual_local_diagonal_sector",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    route_checks = [
        {
            "id": "p22_residual_diagonal_gap_was_still_missing",
            "actual": prior_gap in p22["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P22 still had the residual local diagonal gap as a missing object",
        },
        {
            "id": "r16_declared_control_pullback_packet_present",
            "actual": r16["result"]["explicit_declared_control_pullback_of_residual_local_diagonal_sector_present"],
            "expected": True,
            "meaning": "R16 adds the explicit declared control pullback of the residual local diagonal sector",
        },
        {
            "id": "r16_pair1_declared_block_present",
            "actual": r16["result"]["pair1_declared_residual_diagonal_block_present"],
            "expected": True,
            "meaning": "R16 makes the pair1 declared residual block explicit",
        },
        {
            "id": "r16_zero_or_host_side_cancellation_still_absent",
            "actual": r16["result"]["zero_or_host_side_cancellation_witness_present"],
            "expected": False,
            "meaning": "R16 does not claim zero/cancellation of the residual declared control pullback",
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
        "zero_or_host_side_cancellation_witness_present": False,
        "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
    }

    report = {
        "stage": "P23",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_matching_witness_after_R16",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_WITNESS_ROUTE_AFTER_R16_RESIDUAL_DIAGONAL_PULLBACK_PACKET",
        "reason": "the repo now contains the explicit declared control pullback of the residual local diagonal sector, but it still lacks a witness that this declared pullback vanishes or matches a host-side correction, and the QW-2191 canonicalization boundary remains open",
        "lane": "existing_kernel_feedback_host_matching_witness_route_after_R16",
        "route_under_test": [
            "partial_host_to_canonical_block_overlap_packet",
            "kernel_coefficient_specialization_witness",
            "host_scalar_floor_embedding_packet",
            "declared_control_pullback_of_residual_local_diagonal_sector",
            "zero_or_host_side_cancellation_witness_for_that_pullback",
            "full_physical_uniqueness_or_selector_relevant_canonicalization",
            "host_to_concrete_Psi_block_identification",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R16_explicit_residual_local_diagonal_declared_control_pullback_packet",
            "R15_explicit_host_scalar_floor_embedding_packet",
            "R14_explicit_frozen_kernel_specialization_packet",
            "R13_partial_host_to_canonical_block_overlap_packet",
        ],
        "decomposition_of_P22_residual_diagonal_gap": {
            "from_P22": prior_gap,
            "discharged_by_R16_as": "explicit_declared_control_pullback_of_the_residual_local_diagonal_sector",
            "remaining_matching_gap": "explicit_zero_or_host_side_cancellation_witness_for_the_declared_control_pullback_of_the_residual_local_diagonal_sector",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R16_B1": r16["frontier_text"],
            "QW2191_required_next_step": q2191["required_next_step"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_SELECTOR_RELEVANT_CANONICALIZATION_AND_ADD_A_ZERO_OR_HOST_SIDE_CANCELLATION_WITNESS_FOR_THE_RESIDUAL_DECLARED_PULLBACK_OR_KEEP_THE_HOST_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P23",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P22_residual_diagonal_gap": report["decomposition_of_P22_residual_diagonal_gap"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
