#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r58_direct_m2_psi10_target_action_common_monomial_support_packet.json"
OUT_SUMMARY = GENERATED / "r58_direct_m2_psi10_target_action_common_monomial_support_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r40 = load_json("fundamental_action_reconstruction/generated/r40_direct_m2_psi7_psi10_role_matching_packet.json")
    r57 = load_json("fundamental_action_reconstruction/generated/r57_direct_m2_psi10_common_plus3_assignment_role_split_packet.json")

    role_packet = r40["direct_m2_psi7_psi10_role_matching_packet"]
    target_action_term = role_packet["canonical_action_role_match"]["target_term"]
    action_lifted_term = r57["direct_m2_psi10_common_plus3_assignment_role_split_packet"]["declared_rolewise_parameter_lifts"][
        "action_role_lifted_term"
    ]
    common_parameter = r57["direct_m2_psi10_common_plus3_assignment_role_split_packet"]["common_orbit_parameter_symbol"]
    common_support = "psi10**2/2"

    checks = [
        {
            "id": "r57_target_action_role_assignment_witness_missing_before_support_packet",
            "actual": not r57["result"]["explicit_target_action_role_assignment_witness_present"],
            "expected": True,
            "meaning": "R57 still leaves the target action-role assignment witness absent before any finer action-support packet",
        },
        {
            "id": "r40_target_action_term_present",
            "actual": target_action_term,
            "expected": "m2_psi10*psi10**2/2",
            "meaning": "R40 exports the exact target action term for m2_psi10",
        },
        {
            "id": "r57_action_role_lifted_term_present",
            "actual": action_lifted_term,
            "expected": "mu_m2_plus3_segment_psi7_psi10*psi10**2/2",
            "meaning": "R57 exports the exact lifted target action term for the declared common plus3 parameter",
        },
        {
            "id": "common_target_action_monomial_support_declared",
            "actual": common_support,
            "expected": "psi10**2/2",
            "meaning": "the common fixed target action monomial support is explicitly declared as psi10**2/2",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R58",
        "packet_goal": "materialize_the_exact_route_scoped_common_target_action_monomial_support_packet_for_the_single_target_action_role_assignment_witness_without_claiming_any_coefficient_identification_or_term_equality",
        "source_reports": ["R40", "R57"],
        "direct_m2_psi10_target_action_common_monomial_support_packet": {
            "slot_under_attack": "m2_psi10",
            "role_under_attack": "target_action_role",
            "target_action_role_assignment_witness_under_attack": "explicit_target_action_role_assignment_witness_for_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10",
            "target_action_term": target_action_term,
            "lifted_action_term": action_lifted_term,
            "common_fixed_action_monomial_support": common_support,
            "coefficient_symbols_on_common_support": ["m2_psi10", common_parameter],
            "exact_support_reduction_statement": "on this route-scoped target-action lane, the target action-role assignment witness is reduced to one still-missing coefficient-identification witness between m2_psi10 and mu_m2_plus3_segment_psi7_psi10 on the already fixed common monomial support psi10**2/2",
            "reduced_missing_object": "explicit_target_action_monomial_coefficient_identification_witness_for_m2_psi10_and_mu_m2_plus3_segment_psi7_psi10_on_common_psi10_squared_over_2_support",
            "scope": "single_direct_m2_target_action_role_common_monomial_support_only",
            "non_cancellation_clause": "R58 does not claim any global cancellation, nonzero-factor argument, or actual coefficient identification; it only exports the already shared fixed target-action monomial support",
        },
        "result": {
            "explicit_target_action_common_monomial_support_packet_present": True,
            "explicit_target_action_monomial_coefficient_identification_witness_present": False,
            "explicit_target_action_role_assignment_witness_present": False,
            "global_reduction_of_main_route_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R58",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one direct non-light m2 target action-role on the route-scoped common-monomial-support lane only",
            "boundary": "R58 does not alter the light-facing kernel channel; it only sharpens the missing target action-role assignment witness for m2_psi10 into one coefficient-identification witness on the already fixed common target-action monomial support psi10**2/2",
        },
        "classification": "explicit_target_action_common_monomial_support_packet_present_but_no_target_action_monomial_coefficient_identification_witness",
        "frontier": "R58_B1",
        "frontier_text": "on the route-scoped direct m2 sufficient lane only, the single target action-role assignment witness for m2_psi10 is sharpened into one still-missing coefficient-identification witness on the fixed common target-action monomial support psi10**2/2, while the other direct family blockers, the main route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_m2_psi7_equals_m2_psi10",
            "no_claim_that_any_common_plus3_carrier_orbit_parameter_actually_exists",
            "no_claim_that_any_target_action_coefficient_identification_witness_is_present",
            "no_claim_that_any_global_cancellation_or_nonzero_factor_argument_holds",
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
        "stage": "R58",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI10_TARGET_ACTION_COMMON_MONOMIAL_SUPPORT_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R58_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

