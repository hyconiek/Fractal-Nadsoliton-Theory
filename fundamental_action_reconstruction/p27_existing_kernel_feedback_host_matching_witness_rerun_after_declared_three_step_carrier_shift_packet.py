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
    / "p27_existing_kernel_feedback_host_matching_witness_rerun_after_declared_three_step_carrier_shift_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p27_existing_kernel_feedback_host_matching_witness_rerun_after_declared_three_step_carrier_shift_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p26 = load_json(
        "fundamental_action_reconstruction/generated/p26_existing_kernel_feedback_host_matching_witness_rerun_after_pair1_c1s1_balance_reduction_packet.json"
    )
    r20 = load_json(
        "fundamental_action_reconstruction/generated/r20_declared_three_step_carrier_shift_packet_for_pair1_c1s1_balance_route.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    prior_gap = "explicit_balance_witness_for_the_declared_pair1_residual_c1s1_equation"
    remaining_missing = [
        "explicit_declared_plus3_shift_equivariance_witness_for_the_pair1_c1s1_support_sum",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    route_checks = [
        {
            "id": "p26_c1s1_balance_witness_gap_was_still_missing",
            "actual": prior_gap in p26["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P26 still had the c1s1 balance witness as a missing object",
        },
        {
            "id": "r20_declared_plus3_shift_packet_present",
            "actual": r20["result"]["explicit_declared_plus3_carrier_shift_packet_present"],
            "expected": True,
            "meaning": "R20 adds the explicit declared +3 carrier shift packet",
        },
        {
            "id": "r20_shift_equivariance_witness_still_absent",
            "actual": r20["result"]["explicit_declared_plus3_shift_equivariance_witness_for_pair1_c1s1_support_present"],
            "expected": False,
            "meaning": "R20 does not claim that the c1s1 support sum is invariant under the declared +3 shift",
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
        "pair1_c1s1_balance_reduction_present": True,
        "declared_plus3_carrier_shift_packet_present": True,
        "pair1_c1s1_shift_equivariance_witness_present": False,
        "pair1_c1c1_zero_witness_present": False,
        "pair1_s1s1_zero_witness_present": False,
        "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
    }

    report = {
        "stage": "P27",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_matching_witness_after_R20",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_WITNESS_ROUTE_AFTER_R20_DECLARED_THREE_STEP_CARRIER_SHIFT_PACKET",
        "reason": "the c1s1 balance witness is now reduced to one declared +3 shift-equivariance witness on its support sum, but the repo still proves neither that equivariance nor the c1c1 and s1s1 zero equations, and still leaves the QW-2191 canonicalization boundary open",
        "lane": "existing_kernel_feedback_host_matching_witness_route_after_R20",
        "route_under_test": [
            "partial_host_to_canonical_block_overlap_packet",
            "kernel_coefficient_specialization_witness",
            "host_scalar_floor_embedding_packet",
            "declared_control_pullback_of_residual_local_diagonal_sector",
            "host_side_residual_diagonal_correction_branch",
            "pair1_coefficient_class_reduction",
            "pair1_c1s1_balance_reduction",
            "declared_plus3_carrier_shift_packet",
            "explicit_declared_plus3_shift_equivariance_witness_for_the_pair1_c1s1_support_sum",
            "explicit_zero_witness_for_the_declared_pair1_c1c1_equation",
            "explicit_zero_witness_for_the_declared_pair1_s1s1_equation",
            "full_physical_uniqueness_or_selector_relevant_canonicalization",
            "host_to_concrete_Psi_block_identification",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R20_explicit_declared_plus3_carrier_shift_packet",
            "R19_explicit_declared_pair1_c1s1_balance_reduction_packet",
            "R18_explicit_pair1_coefficient_class_reduction_packet",
            "R14_explicit_frozen_kernel_specialization_packet",
        ],
        "decomposition_of_P26_c1s1_balance_gap": {
            "from_P26": prior_gap,
            "reduced_by_R20_into": "explicit_declared_plus3_shift_equivariance_witness_for_the_pair1_c1s1_support_sum",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R20_B1": r20["frontier_text"],
            "QW2191_required_next_step": q2191["required_next_step"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_THE_DECLARED_PLUS3_SHIFT_EQUIVARIANCE_FOR_THE_PAIR1_C1S1_SUPPORT_SUM_AND_THE_C1C1_S1S1_ZERO_EQUATIONS_AND_DISCHARGE_SELECTOR_RELEVANT_CANONICALIZATION_OR_KEEP_THE_HOST_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P27",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P26_c1s1_balance_gap": report["decomposition_of_P26_c1s1_balance_gap"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
