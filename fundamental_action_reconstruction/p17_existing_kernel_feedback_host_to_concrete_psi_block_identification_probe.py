#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p17_existing_kernel_feedback_host_to_concrete_psi_block_identification_probe.json"
OUT_SUMMARY = GENERATED / "p17_existing_kernel_feedback_host_to_concrete_psi_block_identification_probe_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p16 = load_json(
        "fundamental_action_reconstruction/generated/p16_existing_kernel_feedback_legacy_chart_reduced_operator_export_probe.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")
    r14 = load_json(
        "fundamental_action_reconstruction/generated/r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route_summary.json"
    )
    r15 = load_json(
        "fundamental_action_reconstruction/generated/r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route_summary.json"
    )
    r16 = load_json(
        "fundamental_action_reconstruction/generated/r16_explicit_residual_local_diagonal_declared_control_pullback_packet_for_host_matching_route_summary.json"
    )
    r17 = load_json(
        "fundamental_action_reconstruction/generated/r17_explicit_host_side_residual_diagonal_correction_absence_packet_for_host_matching_route_summary.json"
    )
    r18 = load_json(
        "fundamental_action_reconstruction/generated/r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route_summary.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )

    p16_missing = p16["remaining_missing_upstream_objects"]
    remaining_missing = [
        "explicit_zero_or_host_side_cancellation_witness_for_the_declared_control_pullback_of_the_residual_local_diagonal_sector",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    route_checks = [
        {
            "id": "p16_missing_is_now_the_residual_diagonal_cancellation_and_qw2191_canonicalization_pair",
            "pass": all(obj in p16_missing for obj in remaining_missing),
            "expected": True,
            "actual": all(obj in p16_missing for obj in remaining_missing),
            "meaning": "P16 now reduces the host-to-canonical matching gap to the residual diagonal cancellation witness plus the QW-2191 canonicalization boundary",
        },
        {
            "id": "r14_kernel_channel_specialization_packet_present",
            "pass": r14["status"] == "PASS_PARTIAL_EXPLICIT_FROZEN_KERNEL_SPECIALIZATION_PACKET_READY",
            "expected": "PASS_PARTIAL_EXPLICIT_FROZEN_KERNEL_SPECIALIZATION_PACKET_READY",
            "actual": r14["status"],
            "meaning": "the shared kernel/light-facing channel specialization witness is exported (R14)",
        },
        {
            "id": "r15_scalar_floor_embedding_packet_present",
            "pass": r15["status"] == "PASS_PARTIAL_EXPLICIT_HOST_SCALAR_FLOOR_EMBEDDING_PACKET_READY",
            "expected": "PASS_PARTIAL_EXPLICIT_HOST_SCALAR_FLOOR_EMBEDDING_PACKET_READY",
            "actual": r15["status"],
            "meaning": "the host scalar-floor embedding into the canonical diagonal sector is exported (R15)",
        },
        {
            "id": "r16_residual_diagonal_declared_pullback_packet_present",
            "pass": r16["status"]
            == "PASS_PARTIAL_EXPLICIT_RESIDUAL_DIAGONAL_DECLARED_CONTROL_PULLBACK_PACKET_READY",
            "expected": "PASS_PARTIAL_EXPLICIT_RESIDUAL_DIAGONAL_DECLARED_CONTROL_PULLBACK_PACKET_READY",
            "actual": r16["status"],
            "meaning": "the declared control pullback of the residual local diagonal sector is exported (R16)",
        },
        {
            "id": "r17_host_side_residual_diagonal_correction_absence_packet_present",
            "pass": r17["status"]
            == "PASS_PARTIAL_EXPLICIT_HOST_SIDE_RESIDUAL_DIAGONAL_CORRECTION_ABSENCE_PACKET_READY",
            "expected": "PASS_PARTIAL_EXPLICIT_HOST_SIDE_RESIDUAL_DIAGONAL_CORRECTION_ABSENCE_PACKET_READY",
            "actual": r17["status"],
            "meaning": "the host-side residual diagonal correction branch is closed as absent (R17)",
        },
        {
            "id": "r18_pair1_residual_zero_system_present",
            "pass": r18["status"]
            == "PASS_PARTIAL_PAIR1_RESIDUAL_DECLARED_PULLBACK_COEFFICIENT_CLASS_REDUCTION_PACKET_READY",
            "expected": "PASS_PARTIAL_PAIR1_RESIDUAL_DECLARED_PULLBACK_COEFFICIENT_CLASS_REDUCTION_PACKET_READY",
            "actual": r18["status"],
            "meaning": "the declared pair1 residual block is reduced to an explicit finite zero system (R18), but no zero witness is exported",
        },
        {
            "id": "qw2191_full_physical_uniqueness_still_open",
            "pass": q2191["flags"]["full_physical_uniqueness_closed"] is False,
            "expected": False,
            "actual": q2191["flags"]["full_physical_uniqueness_closed"],
            "meaning": "QW-2191 still blocks full physical uniqueness / selector-relevant canonicalization",
        },
        {
            "id": "c10_host_to_concrete_psi_block_identification_not_shown",
            "pass": c10["result"]["host_to_concrete_psi_block_identification"] == "not_shown",
            "expected": "not_shown",
            "actual": c10["result"]["host_to_concrete_psi_block_identification"],
            "meaning": "the host-to-concrete Psi-block identification is still absent at audit level",
        },
    ]

    route_state = {
        "partial_host_to_canonical_block_overlap_present": True,
        "kernel_channel_specialization_witness_present": True,
        "host_scalar_floor_embedding_present": True,
        "declared_control_pullback_of_residual_local_diagonal_sector_present": True,
        "explicit_pair1_residual_zero_system_present": True,
        "zero_or_host_side_cancellation_witness_present": False,
        "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
    }

    report = {
        "stage": "P17",
        "goal": "compute_or_fail_existing_kernel_feedback_host_to_concrete_Psi_sector_block_identification",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_ROUTE",
        "reason": "the repo now exports a partial host-to-canonical overlap packet (kernel specialization + scalar floor embedding) and the declared control pullback of the residual local diagonal sector reduced to an explicit finite pair1 zero system (R14–R18), but it still exports neither a zero/cancellation witness for that declared residual pullback nor selector-relevant physical canonicalization within the QW-2191 O(2) family; therefore host-to-concrete Psi-block identification remains noncomputable in strict scope",
        "lane": "existing_kernel_feedback_host_to_concrete_Psi_block_identification_route_after_P16",
        "route_under_test": [
            "kernel_channel_specialization_witness",
            "host_scalar_floor_embedding_packet",
            "declared_control_pullback_of_residual_local_diagonal_sector",
            "zero_or_host_side_cancellation_witness_for_that_pullback",
            "full_physical_uniqueness_or_selector_relevant_canonicalization",
            "host_to_concrete_Psi_block_identification",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R18_pair1_residual_declared_pullback_coefficient_class_reduction_packet",
            "R17_host_side_residual_diagonal_correction_absence_packet",
            "R16_explicit_residual_local_diagonal_declared_control_pullback_packet",
            "R15_explicit_host_scalar_floor_embedding_packet",
            "R14_explicit_frozen_kernel_specialization_packet",
        ],
        "decomposition_of_P16_missing_object": {
            "from_P16": "host_to_canonical_matching_gap_reduced_to_residual_diagonal_cancellation_and_qw2191_canonicalization",
            "into_current_blockers": remaining_missing,
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "C10_B1": c10["residual_blockers"]["C10_B1"],
            "R18_B1": r18["result"],
            "QW2191_required_next_step": q2191["required_next_step"],
        },
        "required_next_step": "EITHER_PROVE_SELECTOR_RELEVANT_CANONICALIZATION_AND_ADD_A_ZERO_OR_HOST_SIDE_CANCELLATION_WITNESS_FOR_THE_RESIDUAL_DECLARED_PULLBACK_OR_KEEP_THE_HOST_IDENTIFICATION_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P17",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P16_missing_object": report["decomposition_of_P16_missing_object"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
