#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r62_direct_m2_psi2_common_plus3_assignment_role_split_packet.json"
OUT_SUMMARY = GENERATED / "r62_direct_m2_psi2_common_plus3_assignment_role_split_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r49 = load_json("fundamental_action_reconstruction/generated/r49_direct_m2_psi2_psi5_role_matching_packet.json")
    r53 = load_json(
        "fundamental_action_reconstruction/generated/r53_direct_m2_psi2_psi5_common_plus3_assignment_slot_split_packet.json"
    )

    role_packet = r49["direct_m2_psi2_psi5_role_matching_packet"]
    action_source_term = role_packet["canonical_action_role_match"]["source_term"]
    eom_source_term = role_packet["canonical_eom_role_match"]["source_term"]
    common_parameter = r53["direct_m2_psi2_psi5_common_plus3_assignment_slot_split_packet"][
        "common_orbit_parameter_symbol"
    ]
    source_assignment_witness = "explicit_assignment_witness_of_m2_psi2_to_mu_m2_plus3_segment_psi2_psi5"
    action_lifted_term = f"{common_parameter}*psi2**2/2"
    eom_lifted_term = f"{common_parameter}*psi2(x)"

    checks = [
        {
            "id": "r53_source_slot_assignment_witness_missing_before_split",
            "actual": not r53["result"]["explicit_assignment_witness_of_m2_psi2_to_mu_m2_plus3_segment_psi2_psi5_present"],
            "expected": True,
            "meaning": "R53 still leaves the source-slot assignment witness absent before any finer role split",
        },
        {
            "id": "r49_source_action_role_present",
            "actual": action_source_term,
            "expected": "m2_psi2*psi2**2/2",
            "meaning": "R49 exports the exact source action role for m2_psi2",
        },
        {
            "id": "r49_source_eom_role_present",
            "actual": eom_source_term,
            "expected": "m2_psi2*psi2(x)",
            "meaning": "R49 exports the exact source eom role for m2_psi2",
        },
        {
            "id": "r53_common_parameter_symbol_present",
            "actual": common_parameter,
            "expected": "mu_m2_plus3_segment_psi2_psi5",
            "meaning": "R53 preserves the declared common plus3 carrier-segment parameter symbol on the one-pair route",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R62",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-R53",
        "packet_goal": "materialize_the_exact_route_scoped_role_split_of_the_single_source_slot_assignment_witness_for_m2_psi2_into_action_and_eom_role_assignment_witnesses_without_claiming_either_assignment",
        "source_reports": ["R49", "R53"],
        "direct_m2_psi2_common_plus3_assignment_role_split_packet": {
            "slot_under_attack": "m2_psi2",
            "combined_source_slot_assignment_witness_under_attack": source_assignment_witness,
            "common_orbit_parameter_symbol": common_parameter,
            "source_role_support": {
                "canonical_action_source_term": action_source_term,
                "canonical_eom_source_term": eom_source_term,
            },
            "declared_rolewise_parameter_lifts": {
                "action_role_lifted_term": action_lifted_term,
                "eom_role_lifted_term": eom_lifted_term,
            },
            "exact_role_split": [
                "explicit_source_action_role_assignment_witness_for_m2_psi2_to_mu_m2_plus3_segment_psi2_psi5",
                "explicit_source_eom_role_assignment_witness_for_m2_psi2_to_mu_m2_plus3_segment_psi2_psi5",
            ],
            "exact_split_statement": "on this route-scoped one-pair lane, the source-slot assignment witness for m2_psi2 can only appear through both the source action-role assignment witness and the source eom-role assignment witness to the same declared common plus3 carrier-segment parameter",
            "scope": "single_direct_m2_source_slot_m2_psi2_common_plus3_assignment_role_split_only",
            "non_promotion_clause": "R62 does not claim that either role-specific assignment witness is present in the current repo, nor that the declared common plus3 carrier-segment parameter exists outside the already declared sufficient route",
        },
        "result": {
            "explicit_source_slot_assignment_role_split_packet_present": True,
            "explicit_source_action_role_assignment_witness_present": False,
            "explicit_source_eom_role_assignment_witness_present": False,
            "explicit_source_slot_assignment_witness_present": False,
            "global_reduction_of_main_route_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R62",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one direct non-light m2 source slot on the route-scoped assignment-role-splitting lane only",
            "boundary": "R62 does not alter the light-facing kernel channel; it only sharpens the missing source-slot assignment witness for m2_psi2 into source action-role and source eom-role assignment witnesses",
        },
        "classification": "explicit_source_slot_assignment_role_split_packet_present_but_both_role_specific_assignment_witnesses_absent",
        "frontier": "R62_B1",
        "frontier_text": "on the route-scoped direct m2 sufficient lane only, the single source-slot assignment witness for m2_psi2 is sharpened into two still-missing role-specific assignment witnesses on the source action and source eom terms, while the other direct family blockers, the main route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_m2_psi2_equals_m2_psi5",
            "no_claim_that_any_common_plus3_carrier_orbit_parameter_actually_exists",
            "no_claim_that_either_source_role_assignment_witness_is_present",
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
        "stage": "R62",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI2_COMMON_PLUS3_ASSIGNMENT_ROLE_SPLIT_PACKET_READY",
        "lane": artifact["lane"],
        "result": artifact["classification"],
        "frontier": ["R62_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

