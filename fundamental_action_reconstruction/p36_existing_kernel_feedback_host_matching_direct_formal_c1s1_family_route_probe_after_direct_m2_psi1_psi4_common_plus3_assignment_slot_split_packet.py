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
    / "p36_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_psi4_common_plus3_assignment_slot_split_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p36_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_psi4_common_plus3_assignment_slot_split_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p35 = load_json(
        "fundamental_action_reconstruction/generated/p35_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_psi4_common_plus3_orbit_parameter_packet.json"
    )
    r29 = load_json(
        "fundamental_action_reconstruction/generated/r29_direct_m2_psi1_psi4_common_plus3_assignment_slot_split_packet.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    prior_gap = "explicit_assignment_witness_of_m2_psi1_and_m2_psi4_to_one_common_plus3_carrier_segment_parameter"
    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_assignment_witness_of_m2_psi1_to_mu_m2_plus3_segment_psi1_psi4",
        "explicit_assignment_witness_of_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
        "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10",
        "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
        "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    route_checks = [
        {
            "id": "p35_combined_assignment_gap_was_still_missing",
            "actual": prior_gap in p35["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P35 still had the combined assignment witness for m2_psi1 / m2_psi4 as a missing object",
        },
        {
            "id": "r29_assignment_slot_split_packet_present",
            "actual": r29["result"]["explicit_common_plus3_assignment_slot_split_packet_present"],
            "expected": True,
            "meaning": "R29 adds the exact slotwise split of the one-pair combined assignment witness",
        },
        {
            "id": "r29_source_slot_assignment_witness_still_absent",
            "actual": r29["result"]["explicit_assignment_witness_of_m2_psi1_to_mu_m2_plus3_segment_psi1_psi4_present"],
            "expected": False,
            "meaning": "R29 does not claim that the source-slot assignment witness holds",
        },
        {
            "id": "r29_target_slot_assignment_witness_still_absent",
            "actual": r29["result"]["explicit_assignment_witness_of_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4_present"],
            "expected": False,
            "meaning": "R29 does not claim that the target-slot assignment witness holds",
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

    report = {
        "stage": "P36",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_matching_on_the_direct_formal_pair1_c1s1_family_route_after_R29_direct_m2_psi1_psi4_common_plus3_assignment_slot_split_packet",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_AFTER_R29_DIRECT_M2_PSI1_PSI4_COMMON_PLUS3_ASSIGNMENT_SLOT_SPLIT_PACKET",
        "reason": "on the route-scoped direct m2 sufficient lane only, the single combined assignment witness for m2_psi1 = m2_psi4 is reduced to two still-missing slotwise assignment witnesses to the declared common plus3 carrier-segment parameter, but both slotwise witnesses are absent; the direct g4, g6, and gY family zero witnesses remain absent, the other direct m2 pairwise witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "lane": "current_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_after_R29_only",
        "route_under_test": [
            "main_host_matching_route_after_R21",
            "direct_formal_pair1_c1s1_coefficient_family_route_packet",
            "direct_g4_family_zero_witness",
            "direct_g6_family_zero_witness",
            "direct_gY_family_zero_witness",
            "direct_m2_pairwise_matching_sufficient_route_packet",
            "common_plus3_carrier_segment_parameter_sufficient_route_for_m2_psi1_and_m2_psi4",
            "explicit_assignment_witness_of_m2_psi1_to_mu_m2_plus3_segment_psi1_psi4",
            "explicit_assignment_witness_of_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
            "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10",
            "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
            "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
            "explicit_zero_witness_for_the_declared_pair1_c1c1_equation",
            "explicit_zero_witness_for_the_declared_pair1_s1s1_equation",
            "full_physical_uniqueness_or_selector_relevant_canonicalization",
            "host_to_concrete_Psi_block_identification",
        ],
        "route_checks": route_checks,
        "route_state": {
            "main_route_after_R21_still_negative": True,
            "direct_formal_c1s1_family_route_packet_present": True,
            "direct_g4_family_zero_witness_present": False,
            "direct_g6_family_zero_witness_present": False,
            "direct_gY_family_zero_witness_present": False,
            "common_plus3_carrier_segment_parameter_sufficient_route_present": True,
            "explicit_assignment_witness_of_m2_psi1_to_mu_m2_plus3_segment_psi1_psi4_present": False,
            "explicit_assignment_witness_of_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4_present": False,
            "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10_present": False,
            "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5_present": False,
            "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11_present": False,
            "pair1_c1c1_zero_witness_present": False,
            "pair1_s1s1_zero_witness_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "supporting_present_but_insufficient_objects": [
            "R29_direct_m2_psi1_psi4_common_plus3_assignment_slot_split_packet",
            "R28_direct_m2_psi1_psi4_common_plus3_orbit_parameter_sufficient_route_packet",
            "R27_direct_m2_psi1_psi4_declared_formal_slot_separation_packet",
            "R14_explicit_frozen_kernel_specialization_packet",
        ],
        "decomposition_of_P35_single_pair_gap": {
            "from_P35": prior_gap,
            "reduced_by_R29_into": [
                "explicit_assignment_witness_of_m2_psi1_to_mu_m2_plus3_segment_psi1_psi4",
                "explicit_assignment_witness_of_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
            ],
            "scope": "single_direct_m2_pair_m2_psi1_and_m2_psi4_common_plus3_carrier_segment_parameter_slot_split_only",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R29_B1": r29["frontier_text"],
            "QW2191_required_next_step": q2191["required_next_step"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_ONE_OF_THE_TWO_SLOTWISE_ASSIGNMENT_WITNESSES_FOR_M2_PSI1_OR_M2_PSI4_TO_THE_DECLARED_COMMON_PLUS3_CARRIER_SEGMENT_PARAMETER_ON_THIS_ROUTE_OR_ATTACK_ONE_OF_THE_OTHER_DIRECT_M2_PAIRWISE_WITNESSES_OR_ATTACK_ONE_OF_THE_DIRECT_G4_G6_GY_ZERO_WITNESSES_AND_PROVE_THE_PAIR1_C1C1_AND_S1S1_ZERO_EQUATIONS_AND_DISCHARGE_SELECTOR_RELEVANT_CANONICALIZATION_OR_KEEP_BOTH_THE_MAIN_AND_DIRECT_ROUTES_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P36",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P35_single_pair_gap": report["decomposition_of_P35_single_pair_gap"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
