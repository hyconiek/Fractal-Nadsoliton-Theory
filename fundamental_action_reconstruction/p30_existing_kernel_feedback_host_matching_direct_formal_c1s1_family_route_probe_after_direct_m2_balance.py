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
    / "p30_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_balance.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p30_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_balance_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p29 = load_json(
        "fundamental_action_reconstruction/generated/p29_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe.json"
    )
    r23 = load_json(
        "fundamental_action_reconstruction/generated/r23_direct_mass_like_m2_family_balance_reduction_packet.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    prior_gap = "explicit_zero_witness_for_direct_mass_like_m2_family_c1s1_shift_defect"
    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_balance_witness_for_direct_mass_like_m2_family_c1s1_shift_defect",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    route_checks = [
        {
            "id": "p29_direct_m2_zero_witness_was_still_missing",
            "actual": prior_gap in p29["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P29 still had the direct mass-like m2 family zero witness as a missing object",
        },
        {
            "id": "r23_direct_m2_balance_reduction_present",
            "actual": r23["result"]["explicit_direct_mass_like_m2_family_balance_reduction_present"],
            "expected": True,
            "meaning": "R23 adds the exact direct mass-like m2 family balance reduction",
        },
        {
            "id": "r23_direct_m2_balance_witness_still_absent",
            "actual": r23["result"]["explicit_direct_mass_like_m2_family_balance_witness_present"],
            "expected": False,
            "meaning": "R23 does not claim that the direct mass-like m2 family balance holds",
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
        "stage": "P30",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_matching_on_the_direct_formal_pair1_c1s1_family_route_after_R23_direct_m2_balance_reduction",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_AFTER_R23_DIRECT_M2_BALANCE_PACKET",
        "reason": "on the direct formal coefficient-family route only, the direct mass-like m2 family zero witness is now reduced to one exact balance witness, but that balance is not proved; the direct g4, g6, and gY family zero witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "lane": "current_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_after_R23_only",
        "route_under_test": [
            "main_host_matching_route_after_R21",
            "direct_formal_pair1_c1s1_coefficient_family_route_packet",
            "direct_g4_family_zero_witness",
            "direct_g6_family_zero_witness",
            "direct_gY_family_zero_witness",
            "direct_mass_like_m2_family_balance_reduction",
            "direct_mass_like_m2_family_balance_witness",
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
            "direct_mass_like_m2_family_balance_reduction_present": True,
            "direct_mass_like_m2_family_balance_witness_present": False,
            "pair1_c1c1_zero_witness_present": False,
            "pair1_s1s1_zero_witness_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "supporting_present_but_insufficient_objects": [
            "R23_direct_mass_like_m2_family_balance_reduction_packet",
            "R22_direct_formal_c1s1_family_route_packet",
            "R21_explicit_pair1_c1s1_shift_defect_polynomial_packet",
            "R14_explicit_frozen_kernel_specialization_packet",
        ],
        "decomposition_of_P29_direct_m2_zero_gap": {
            "from_P29": prior_gap,
            "reduced_by_R23_into": "explicit_balance_witness_for_direct_mass_like_m2_family_c1s1_shift_defect",
            "scope": "direct_mass_like_m2_family_pair1_c1s1_route_only",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R23_B1": r23["frontier_text"],
            "QW2191_required_next_step": q2191["required_next_step"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_THE_DIRECT_M2_BALANCE_AND_THE_DIRECT_G4_G6_GY_ZERO_WITNESSES_AND_PROVE_THE_PAIR1_C1C1_AND_S1S1_ZERO_EQUATIONS_AND_DISCHARGE_SELECTOR_RELEVANT_CANONICALIZATION_OR_KEEP_BOTH_THE_MAIN_AND_DIRECT_ROUTES_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P30",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P29_direct_m2_zero_gap": report["decomposition_of_P29_direct_m2_zero_gap"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
