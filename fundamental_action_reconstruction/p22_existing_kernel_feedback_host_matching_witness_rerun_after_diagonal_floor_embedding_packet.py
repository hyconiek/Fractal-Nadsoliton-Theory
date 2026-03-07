#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p22_existing_kernel_feedback_host_matching_witness_rerun_after_diagonal_floor_embedding_packet.json"
OUT_SUMMARY = GENERATED / "p22_existing_kernel_feedback_host_matching_witness_rerun_after_diagonal_floor_embedding_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p21 = load_json(
        "fundamental_action_reconstruction/generated/p21_existing_kernel_feedback_host_matching_witness_rerun_after_kernel_specialization_packet.json"
    )
    r15 = load_json(
        "fundamental_action_reconstruction/generated/r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    prior_gap = (
        "explicit_diagonal_sector_equality_or_matching_witness_linking_the_host_floor_m0_squared_I_to_the_canonical_local_diagonal_sector_or_to_a_declared_control_pullback_of_it"
    )
    remaining_missing = [
        "explicit_residual_local_diagonal_sector_equality_or_cancellation_witness_reducing_the_canonical_diagonal_sector_to_the_host_floor_m0_squared_I_or_to_a_declared_control_pullback_of_it",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    route_checks = [
        {
            "id": "p21_diagonal_matching_gap_was_still_missing",
            "actual": prior_gap in p21["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P21 still had the diagonal matching witness as a missing object",
        },
        {
            "id": "r15_host_scalar_floor_embedding_packet_present",
            "actual": r15["embedding_result"]["explicit_host_scalar_floor_embedding_packet_present"],
            "expected": True,
            "meaning": "R15 adds the explicit scalar-floor embedding packet",
        },
        {
            "id": "r15_residual_local_diagonal_sector_present",
            "actual": r15["embedding_result"]["residual_local_diagonal_sector_present"],
            "expected": True,
            "meaning": "R15 keeps the residual local diagonal sector explicit",
        },
        {
            "id": "r15_full_diagonal_matching_still_absent",
            "actual": r15["embedding_result"]["full_diagonal_sector_matching_witness_present"],
            "expected": False,
            "meaning": "R15 does not claim full diagonal-sector matching",
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
        "residual_local_diagonal_sector_cancellation_witness_present": False,
        "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
    }

    report = {
        "stage": "P22",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_matching_witness_after_R15",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_WITNESS_ROUTE_AFTER_R15_DIAGONAL_FLOOR_EMBEDDING_PACKET",
        "reason": "the repo now contains an explicit embedding of the host scalar floor into the canonical diagonal sector, but it still lacks a witness that the residual local diagonal sector cancels or equals a declared control pullback, and the QW-2191 canonicalization boundary remains open",
        "lane": "existing_kernel_feedback_host_matching_witness_route_after_R15",
        "route_under_test": [
            "partial_host_to_canonical_block_overlap_packet",
            "kernel_coefficient_specialization_witness",
            "host_scalar_floor_embedding_packet",
            "residual_local_diagonal_sector_cancellation_witness",
            "full_physical_uniqueness_or_selector_relevant_canonicalization",
            "host_to_concrete_Psi_block_identification",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R15_explicit_host_scalar_floor_embedding_packet",
            "R14_explicit_frozen_kernel_specialization_packet",
            "R13_partial_host_to_canonical_block_overlap_packet",
            "R12_explicit_canonical_Psi_block_export",
        ],
        "decomposition_of_P21_diagonal_gap": {
            "from_P21": prior_gap,
            "discharged_by_R15_as": "explicit_host_scalar_floor_embedding_into_the_canonical_diagonal_sector",
            "remaining_matching_gap": "explicit_residual_local_diagonal_sector_equality_or_cancellation_witness_reducing_the_canonical_diagonal_sector_to_the_host_floor_m0_squared_I_or_to_a_declared_control_pullback_of_it",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R15_B1": r15["frontier_text"],
            "QW2191_required_next_step": q2191["required_next_step"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_SELECTOR_RELEVANT_CANONICALIZATION_AND_ADD_A_RESIDUAL_LOCAL_DIAGONAL_CANCELLATION_WITNESS_OR_KEEP_THE_HOST_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P22",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P21_diagonal_gap": report["decomposition_of_P21_diagonal_gap"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
