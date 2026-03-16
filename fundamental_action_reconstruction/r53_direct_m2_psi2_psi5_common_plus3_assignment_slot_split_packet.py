#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r53_direct_m2_psi2_psi5_common_plus3_assignment_slot_split_packet.json"
OUT_SUMMARY = GENERATED / "r53_direct_m2_psi2_psi5_common_plus3_assignment_slot_split_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r52 = load_json(
        "fundamental_action_reconstruction/generated/r52_direct_m2_psi2_psi5_common_plus3_orbit_parameter_sufficient_route_packet.json"
    )

    route = r52["direct_m2_psi2_psi5_common_plus3_orbit_parameter_sufficient_route"]
    assignments = route["sufficient_assignment_route"]
    common_parameter = route["common_orbit_parameter_symbol"]
    source_slot = "m2_psi2"
    target_slot = "m2_psi5"
    source_assignment = f"{source_slot} = {common_parameter}"
    target_assignment = f"{target_slot} = {common_parameter}"

    checks = [
        {
            "id": "r52_common_plus3_route_present",
            "actual": r52["result"]["explicit_common_plus3_carrier_orbit_parameter_sufficient_route_present"],
            "expected": True,
            "meaning": "R52 already exports the one-pair common plus3 carrier-segment parameter sufficient route",
        },
        {
            "id": "r52_source_slot_assignment_present_in_route",
            "actual": source_assignment in assignments,
            "expected": True,
            "meaning": "R52 already exports the source-slot assignment statement inside the combined sufficient route",
        },
        {
            "id": "r52_target_slot_assignment_present_in_route",
            "actual": target_assignment in assignments,
            "expected": True,
            "meaning": "R52 already exports the target-slot assignment statement inside the combined sufficient route",
        },
        {
            "id": "r52_combined_assignment_witness_still_absent",
            "actual": r52["result"]["explicit_assignment_witness_to_one_common_plus3_carrier_orbit_parameter_present"],
            "expected": False,
            "meaning": "R52 does not claim that the combined assignment witness actually holds",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R53",
        "packet_goal": "materialize_the_exact_route_scoped_slotwise_split_of_the_single_direct_m2_pair_combined_assignment_witness_into_two_explicit_slotwise_assignment_witnesses_without_claiming_either_assignment",
        "source_reports": ["R52"],
        "direct_m2_psi2_psi5_common_plus3_assignment_slot_split_packet": {
            "pairwise_target_under_attack": "m2_psi2 = m2_psi5",
            "common_orbit_parameter_symbol": common_parameter,
            "combined_assignment_witness_under_attack": "explicit_assignment_witness_of_m2_psi2_and_m2_psi5_to_one_common_plus3_carrier_segment_parameter",
            "exact_slotwise_split": [
                "explicit_assignment_witness_of_m2_psi2_to_mu_m2_plus3_segment_psi2_psi5",
                "explicit_assignment_witness_of_m2_psi5_to_mu_m2_plus3_segment_psi2_psi5",
            ],
            "slotwise_assignment_statements": [source_assignment, target_assignment],
            "exact_split_statement": "on this route-scoped one-pair sufficient lane, the combined assignment witness is present only if both listed slotwise assignment witnesses are present",
            "scope": "single_direct_m2_pair_m2_psi2_and_m2_psi5_common_plus3_carrier_segment_parameter_slot_split_only",
            "non_promotion_clause": "R53 does not claim that either slotwise assignment witness is present in the current repo, nor that the common plus3 carrier-segment parameter exists outside the already declared sufficient route",
        },
        "result": {
            "explicit_common_plus3_assignment_slot_split_packet_present": True,
            "explicit_assignment_witness_of_m2_psi2_to_mu_m2_plus3_segment_psi2_psi5_present": False,
            "explicit_assignment_witness_of_m2_psi5_to_mu_m2_plus3_segment_psi2_psi5_present": False,
            "explicit_combined_assignment_witness_present": False,
            "global_reduction_of_main_route_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R53",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one direct non-light m2 pair on the route-scoped assignment-splitting lane only",
            "boundary": "R53 does not alter the light-facing kernel channel; it only sharpens the missing combined assignment witness for m2_psi2 / m2_psi5 into two explicit slotwise assignment witnesses",
        },
        "classification": "explicit_common_plus3_assignment_slot_split_packet_present_but_both_slotwise_assignment_witnesses_absent",
        "frontier": "R53_B1",
        "frontier_text": "on the route-scoped direct m2 sufficient lane only, the single combined assignment witness for m2_psi2 / m2_psi5 is sharpened into two still-missing slotwise assignment witnesses to the declared common plus3 carrier-segment parameter, while the other direct family blockers, the main route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_m2_psi2_equals_m2_psi5",
            "no_claim_that_any_common_plus3_carrier_orbit_parameter_actually_exists",
            "no_claim_that_either_slotwise_assignment_witness_is_present",
            "no_claim_that_any_other_direct_m2_pairwise_equality_holds",
            "no_claim_that_the_direct_m2_shift_equivariance_holds",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R53",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI2_PSI5_COMMON_PLUS3_ASSIGNMENT_SLOT_SPLIT_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R53_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

