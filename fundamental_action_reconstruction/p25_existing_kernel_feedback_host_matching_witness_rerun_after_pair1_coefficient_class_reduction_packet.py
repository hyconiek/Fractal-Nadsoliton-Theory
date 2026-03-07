#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p25_existing_kernel_feedback_host_matching_witness_rerun_after_pair1_coefficient_class_reduction_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p25_existing_kernel_feedback_host_matching_witness_rerun_after_pair1_coefficient_class_reduction_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p24 = load_json(
        "fundamental_action_reconstruction/generated/p24_existing_kernel_feedback_host_matching_witness_rerun_after_host_side_residual_absence_packet.json"
    )
    r18 = load_json(
        "fundamental_action_reconstruction/generated/r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    prior_gap = "explicit_zero_witness_for_the_declared_control_pullback_of_the_residual_local_diagonal_sector"
    remaining_missing = [
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1s1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    route_checks = [
        {
            "id": "p24_generic_zero_witness_gap_was_still_missing",
            "actual": prior_gap in p24["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P24 still had the generic zero witness for the declared residual pullback as a missing object",
        },
        {
            "id": "r18_pair1_coefficient_class_reduction_present",
            "actual": r18["result"]["explicit_transport_induced_pair1_coefficient_class_reduction_present"],
            "expected": True,
            "meaning": "R18 adds the exact pair1 coefficient-class reduction",
        },
        {
            "id": "r18_pair1_zero_condition_system_present",
            "actual": r18["result"]["explicit_zero_condition_system_for_declared_pair1_residual_block_present"],
            "expected": True,
            "meaning": "R18 exports the finite exact zero-equation system on pair1",
        },
        {
            "id": "r18_pair1_zero_witness_still_absent",
            "actual": r18["result"]["explicit_zero_witness_for_declared_pair1_residual_block_present"],
            "expected": False,
            "meaning": "R18 does not claim that any pair1 zero equation is satisfied",
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
        "shared_kernel_light_channel_specialized": True,
        "host_scalar_floor_embedding_packet_present": True,
        "declared_control_pullback_of_residual_local_diagonal_sector_present": True,
        "host_side_residual_diagonal_correction_present": False,
        "pair1_coefficient_class_reduction_present": True,
        "pair1_c1c1_zero_witness_present": False,
        "pair1_c1s1_zero_witness_present": False,
        "pair1_s1s1_zero_witness_present": False,
        "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
    }

    report = {
        "stage": "P25",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_matching_witness_after_R18",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_WITNESS_ROUTE_AFTER_R18_PAIR1_COEFFICIENT_CLASS_REDUCTION_PACKET",
        "reason": "the old generic zero-witness blocker is now reduced to three explicit pair1 zero-equation witnesses, but the repo still proves none of them and still leaves the QW-2191 canonicalization boundary open",
        "lane": "existing_kernel_feedback_host_matching_witness_route_after_R18",
        "route_under_test": [
            "partial_host_to_canonical_block_overlap_packet",
            "kernel_coefficient_specialization_witness",
            "host_scalar_floor_embedding_packet",
            "declared_control_pullback_of_residual_local_diagonal_sector",
            "host_side_residual_diagonal_correction_branch",
            "pair1_coefficient_class_reduction",
            "explicit_zero_witness_for_each_independent_pair1_residual_equation",
            "full_physical_uniqueness_or_selector_relevant_canonicalization",
            "host_to_concrete_Psi_block_identification",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R18_explicit_pair1_coefficient_class_reduction_packet",
            "R17_explicit_host_side_residual_diagonal_correction_absence_packet",
            "R16_explicit_residual_local_diagonal_declared_control_pullback_packet",
            "R14_explicit_frozen_kernel_specialization_packet",
        ],
        "decomposition_of_P24_zero_witness_gap": {
            "from_P24": prior_gap,
            "reduced_by_R18_into": remaining_missing[:3],
            "still_open_boundary": remaining_missing[3],
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R18_B1": r18["frontier_text"],
            "QW2191_required_next_step": q2191["required_next_step"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_THE_THREE_INDEPENDENT_DECLARED_PAIR1_ZERO_EQUATIONS_AND_DISCHARGE_SELECTOR_RELEVANT_CANONICALIZATION_OR_KEEP_THE_HOST_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P25",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P24_zero_witness_gap": report["decomposition_of_P24_zero_witness_gap"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
