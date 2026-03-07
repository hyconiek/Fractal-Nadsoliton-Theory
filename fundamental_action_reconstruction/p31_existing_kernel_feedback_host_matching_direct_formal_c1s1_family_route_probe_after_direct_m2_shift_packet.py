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
    / "p31_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_shift_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p31_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_shift_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p30 = load_json(
        "fundamental_action_reconstruction/generated/p30_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_balance.json"
    )
    r24 = load_json(
        "fundamental_action_reconstruction/generated/r24_declared_plus3_shift_packet_for_direct_m2_family_route.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    prior_gap = "explicit_balance_witness_for_direct_mass_like_m2_family_c1s1_shift_defect"
    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_declared_plus3_shift_equivariance_witness_for_direct_mass_like_m2_family_positive_support_sum",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    route_checks = [
        {
            "id": "p30_direct_m2_balance_witness_was_still_missing",
            "actual": prior_gap in p30["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P30 still had the direct m2 balance witness as a missing object",
        },
        {
            "id": "r24_direct_m2_shift_packet_present",
            "actual": r24["result"]["explicit_declared_plus3_shift_packet_for_direct_m2_family_present"],
            "expected": True,
            "meaning": "R24 adds the exact declared +3 shift packet for the direct m2 family route",
        },
        {
            "id": "r24_direct_m2_shift_equivariance_witness_still_absent",
            "actual": r24["result"]["explicit_declared_plus3_shift_equivariance_witness_for_direct_m2_positive_support_present"],
            "expected": False,
            "meaning": "R24 does not claim that the direct m2 shift-equivariance holds",
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
        "stage": "P31",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_matching_on_the_direct_formal_pair1_c1s1_family_route_after_R24_direct_m2_shift_packet",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_AFTER_R24_DIRECT_M2_SHIFT_PACKET",
        "reason": "on the direct formal coefficient-family route only, the direct mass-like m2 family balance witness is now reduced to one declared plus3 shift-equivariance witness, but that shift-equivariance is not proved; the direct g4, g6, and gY family zero witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "lane": "current_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_after_R24_only",
        "route_under_test": [
            "main_host_matching_route_after_R21",
            "direct_formal_pair1_c1s1_coefficient_family_route_packet",
            "direct_g4_family_zero_witness",
            "direct_g6_family_zero_witness",
            "direct_gY_family_zero_witness",
            "declared_plus3_shift_packet_for_direct_m2_family_route",
            "explicit_declared_plus3_shift_equivariance_witness_for_direct_m2_positive_support",
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
            "declared_plus3_shift_packet_for_direct_m2_family_present": True,
            "direct_m2_shift_equivariance_witness_present": False,
            "pair1_c1c1_zero_witness_present": False,
            "pair1_s1s1_zero_witness_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "supporting_present_but_insufficient_objects": [
            "R24_declared_plus3_shift_packet_for_direct_m2_family_route",
            "R23_direct_mass_like_m2_family_balance_reduction_packet",
            "R22_direct_formal_c1s1_family_route_packet",
            "R14_explicit_frozen_kernel_specialization_packet",
        ],
        "decomposition_of_P30_direct_m2_balance_gap": {
            "from_P30": prior_gap,
            "reduced_by_R24_into": "explicit_declared_plus3_shift_equivariance_witness_for_direct_mass_like_m2_family_positive_support_sum",
            "scope": "declared_plus3_shift_on_direct_mass_like_m2_family_positive_support_only",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R24_B1": r24["frontier_text"],
            "QW2191_required_next_step": q2191["required_next_step"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_THE_DIRECT_M2_SHIFT_EQUIVARIANCE_AND_THE_DIRECT_G4_G6_GY_ZERO_WITNESSES_AND_PROVE_THE_PAIR1_C1C1_AND_S1S1_ZERO_EQUATIONS_AND_DISCHARGE_SELECTOR_RELEVANT_CANONICALIZATION_OR_KEEP_BOTH_THE_MAIN_AND_DIRECT_ROUTES_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P31",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P30_direct_m2_balance_gap": report["decomposition_of_P30_direct_m2_balance_gap"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
