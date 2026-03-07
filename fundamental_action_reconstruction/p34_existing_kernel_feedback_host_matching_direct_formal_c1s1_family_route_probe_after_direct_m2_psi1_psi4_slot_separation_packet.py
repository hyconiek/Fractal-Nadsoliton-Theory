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
    / "p34_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_psi4_slot_separation_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p34_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_psi4_slot_separation_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p33 = load_json(
        "fundamental_action_reconstruction/generated/p33_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_psi4_role_matching_packet.json"
    )
    r27 = load_json(
        "fundamental_action_reconstruction/generated/r27_direct_m2_psi1_psi4_declared_formal_slot_separation_packet.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    prior_gap = "explicit_coefficient_identification_witness_for_the_declared_role_matched_direct_m2_pair_m2_psi1_equals_m2_psi4"
    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_common_parameter_source_or_symbol_identification_witness_for_the_declared_formal_m2_slots_m2_psi1_and_m2_psi4",
        "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10",
        "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
        "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    route_checks = [
        {
            "id": "p33_single_coefficient_gap_was_still_missing",
            "actual": prior_gap in p33["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P33 still had the single coefficient-identification witness m2_psi1 = m2_psi4 as a missing object",
        },
        {
            "id": "r27_formal_slot_separation_packet_present",
            "actual": r27["result"]["explicit_declared_formal_slot_separation_packet_for_m2_psi1_and_m2_psi4_present"],
            "expected": True,
            "meaning": "R27 adds the exact one-pair declared formal slot-separation packet",
        },
        {
            "id": "r27_common_source_or_symbol_identification_still_absent",
            "actual": r27["result"]["explicit_common_parameter_source_or_symbol_identification_witness_present"],
            "expected": False,
            "meaning": "R27 does not claim a common parameter-source or symbol-identification witness",
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
        "stage": "P34",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_matching_on_the_direct_formal_pair1_c1s1_family_route_after_R27_direct_m2_psi1_psi4_slot_separation_packet",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_AFTER_R27_DIRECT_M2_PSI1_PSI4_SLOT_SEPARATION_PACKET",
        "reason": "on the route-scoped direct m2 sufficient lane only, the single coefficient-identification gap for m2_psi1 = m2_psi4 is reduced to one still-missing common-parameter-source or symbol-identification witness between two distinct named formal slots, but that witness is absent; the direct g4, g6, and gY family zero witnesses remain absent, the other direct m2 pairwise witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "lane": "current_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_after_R27_only",
        "route_under_test": [
            "main_host_matching_route_after_R21",
            "direct_formal_pair1_c1s1_coefficient_family_route_packet",
            "direct_g4_family_zero_witness",
            "direct_g6_family_zero_witness",
            "direct_gY_family_zero_witness",
            "direct_m2_pairwise_matching_sufficient_route_packet",
            "declared_formal_slot_separation_packet_for_m2_psi1_and_m2_psi4",
            "explicit_common_parameter_source_or_symbol_identification_witness_for_the_declared_formal_m2_slots_m2_psi1_and_m2_psi4",
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
            "declared_formal_slot_separation_packet_for_m2_psi1_and_m2_psi4_present": True,
            "explicit_common_parameter_source_or_symbol_identification_witness_present": False,
            "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10_present": False,
            "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5_present": False,
            "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11_present": False,
            "pair1_c1c1_zero_witness_present": False,
            "pair1_s1s1_zero_witness_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "supporting_present_but_insufficient_objects": [
            "R27_direct_m2_psi1_psi4_declared_formal_slot_separation_packet",
            "R26_direct_m2_psi1_psi4_role_matching_packet",
            "R25_direct_m2_pairwise_matching_sufficient_route_packet",
            "R14_explicit_frozen_kernel_specialization_packet",
        ],
        "decomposition_of_P33_single_pair_gap": {
            "from_P33": prior_gap,
            "reduced_by_R27_into": "explicit_common_parameter_source_or_symbol_identification_witness_for_the_declared_formal_m2_slots_m2_psi1_and_m2_psi4",
            "scope": "single_direct_m2_pair_m2_psi1_and_m2_psi4_declared_formal_slot_separation_only",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R27_B1": r27["frontier_text"],
            "QW2191_required_next_step": q2191["required_next_step"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_THE_SINGLE_DIRECT_M2_COMMON_PARAMETER_SOURCE_OR_SYMBOL_IDENTIFICATION_WITNESS_FOR_M2_PSI1_AND_M2_PSI4_ON_THIS_ROUTE_OR_ATTACK_ONE_OF_THE_OTHER_DIRECT_M2_PAIRWISE_WITNESSES_OR_ATTACK_ONE_OF_THE_DIRECT_G4_G6_GY_ZERO_WITNESSES_AND_PROVE_THE_PAIR1_C1C1_AND_S1S1_ZERO_EQUATIONS_AND_DISCHARGE_SELECTOR_RELEVANT_CANONICALIZATION_OR_KEEP_BOTH_THE_MAIN_AND_DIRECT_ROUTES_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P34",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P33_single_pair_gap": report["decomposition_of_P33_single_pair_gap"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
