#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r57_direct_m2_psi10_common_plus3_assignment_role_split_packet.json"
OUT_SUMMARY = GENERATED / "r57_direct_m2_psi10_common_plus3_assignment_role_split_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r40 = load_json("fundamental_action_reconstruction/generated/r40_direct_m2_psi7_psi10_role_matching_packet.json")
    r43 = load_json(
        "fundamental_action_reconstruction/generated/r43_direct_m2_psi7_psi10_common_plus3_assignment_slot_split_packet.json"
    )

    role_packet = r40["direct_m2_psi7_psi10_role_matching_packet"]
    action_target_term = role_packet["canonical_action_role_match"]["target_term"]
    eom_target_term = role_packet["canonical_eom_role_match"]["target_term"]
    common_parameter = r43["direct_m2_psi7_psi10_common_plus3_assignment_slot_split_packet"]["common_orbit_parameter_symbol"]
    target_assignment_witness = "explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10"
    action_lifted_term = f"{common_parameter}*psi10**2/2"
    eom_lifted_term = f"{common_parameter}*psi10(x)"

    checks = [
        {
            "id": "r43_target_slot_assignment_witness_missing_before_split",
            "actual": not r43["result"]["explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10_present"],
            "expected": True,
            "meaning": "R43 still leaves the target-slot assignment witness absent before any finer role split",
        },
        {
            "id": "r40_target_action_role_present",
            "actual": action_target_term,
            "expected": "m2_psi10*psi10**2/2",
            "meaning": "R40 exports the exact target action role for m2_psi10",
        },
        {
            "id": "r40_target_eom_role_present",
            "actual": eom_target_term,
            "expected": "m2_psi10*psi10(x)",
            "meaning": "R40 exports the exact target eom role for m2_psi10",
        },
        {
            "id": "r43_common_parameter_symbol_present",
            "actual": common_parameter,
            "expected": "mu_m2_plus3_segment_psi7_psi10",
            "meaning": "R43 preserves the declared common plus3 carrier-segment parameter symbol on the one-pair route",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R57",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-AX11-AX12-AX13-R43",
        "packet_goal": "materialize_the_exact_route_scoped_role_split_of_the_single_target_slot_assignment_witness_for_m2_psi10_into_action_and_eom_role_assignment_witnesses_without_claiming_either_assignment",
        "source_reports": ["R40", "R43"],
        "direct_m2_psi10_common_plus3_assignment_role_split_packet": {
            "slot_under_attack": "m2_psi10",
            "combined_target_slot_assignment_witness_under_attack": target_assignment_witness,
            "common_orbit_parameter_symbol": common_parameter,
            "target_role_support": {
                "canonical_action_target_term": action_target_term,
                "canonical_eom_target_term": eom_target_term,
            },
            "declared_rolewise_parameter_lifts": {
                "action_role_lifted_term": action_lifted_term,
                "eom_role_lifted_term": eom_lifted_term,
            },
            "exact_role_split": [
                "explicit_target_action_role_assignment_witness_for_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10",
                "explicit_target_eom_role_assignment_witness_for_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10",
            ],
            "exact_split_statement": "on this route-scoped one-pair lane, the target-slot assignment witness for m2_psi10 can only appear through both the target action-role assignment witness and the target eom-role assignment witness to the same declared common plus3 carrier-segment parameter",
            "scope": "single_direct_m2_target_slot_m2_psi10_common_plus3_assignment_role_split_only",
            "non_promotion_clause": "R57 does not claim that either role-specific assignment witness is present in the current repo, nor that the declared common plus3 carrier-segment parameter exists outside the already declared sufficient route",
        },
        "result": {
            "explicit_target_slot_assignment_role_split_packet_present": True,
            "explicit_target_action_role_assignment_witness_present": False,
            "explicit_target_eom_role_assignment_witness_present": False,
            "explicit_target_slot_assignment_witness_present": False,
            "global_reduction_of_main_route_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R57",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one direct non-light m2 target slot on the route-scoped assignment-role-splitting lane only",
            "boundary": "R57 does not alter the light-facing kernel channel; it only sharpens the missing target-slot assignment witness for m2_psi10 into target action-role and target eom-role assignment witnesses",
        },
        "classification": "explicit_target_slot_assignment_role_split_packet_present_but_both_role_specific_assignment_witnesses_absent",
        "frontier": "R57_B1",
        "frontier_text": "on the route-scoped direct m2 sufficient lane only, the single target-slot assignment witness for m2_psi10 is sharpened into two still-missing role-specific assignment witnesses on the target action and target eom terms, while the other direct family blockers, the main route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_m2_psi7_equals_m2_psi10",
            "no_claim_that_any_common_plus3_carrier_orbit_parameter_actually_exists",
            "no_claim_that_either_target_role_assignment_witness_is_present",
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
        "stage": "R57",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI10_COMMON_PLUS3_ASSIGNMENT_ROLE_SPLIT_PACKET_READY",
        "lane": artifact["lane"],
        "result": artifact["classification"],
        "frontier": ["R57_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

