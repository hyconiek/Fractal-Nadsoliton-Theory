#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r74_direct_m2_psi8_source_action_common_monomial_support_packet.json"
OUT_SUMMARY = GENERATED / "r74_direct_m2_psi8_source_action_common_monomial_support_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r50 = load_json("fundamental_action_reconstruction/generated/r50_direct_m2_psi8_psi11_role_matching_packet.json")
    r72 = load_json("fundamental_action_reconstruction/generated/r72_direct_m2_psi8_common_plus3_assignment_role_split_packet.json")

    role_packet = r50["direct_m2_psi8_psi11_role_matching_packet"]
    action_source_term = role_packet["canonical_action_role_match"]["source_term"]
    action_lifted_term = r72["direct_m2_psi8_common_plus3_assignment_role_split_packet"]["declared_rolewise_parameter_lifts"][
        "action_role_lifted_term"
    ]
    common_parameter = r72["direct_m2_psi8_common_plus3_assignment_role_split_packet"]["common_orbit_parameter_symbol"]
    common_support = "psi8**2/2"

    checks = [
        {
            "id": "r72_source_action_role_assignment_witness_missing_before_support_packet",
            "actual": not r72["result"]["explicit_source_action_role_assignment_witness_present"],
            "expected": True,
            "meaning": "R72 still leaves the source action-role assignment witness absent before any finer action-support packet",
        },
        {
            "id": "r50_source_action_term_present",
            "actual": action_source_term,
            "expected": "m2_psi8*psi8**2/2",
            "meaning": "R50 exports the exact source action term for m2_psi8",
        },
        {
            "id": "r72_action_role_lifted_term_present",
            "actual": action_lifted_term,
            "expected": "mu_m2_plus3_segment_psi8_psi11*psi8**2/2",
            "meaning": "R72 exports the exact lifted source action term for the declared common plus3 parameter",
        },
        {
            "id": "common_action_monomial_support_declared",
            "actual": common_support,
            "expected": "psi8**2/2",
            "meaning": "the common fixed source action monomial support is explicitly declared as psi8**2/2",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R74",
        "packet_goal": "materialize_the_exact_route_scoped_common_source_action_monomial_support_packet_for_the_single_source_action_role_assignment_witness_without_claiming_any_coefficient_identification_or_term_equality",
        "source_reports": ["R50", "R72"],
        "direct_m2_psi8_source_action_common_monomial_support_packet": {
            "slot_under_attack": "m2_psi8",
            "role_under_attack": "source_action_role",
            "source_action_role_assignment_witness_under_attack": "explicit_source_action_role_assignment_witness_for_m2_psi8_to_mu_m2_plus3_segment_psi8_psi11",
            "source_action_term": action_source_term,
            "lifted_action_term": action_lifted_term,
            "common_fixed_action_monomial_support": common_support,
            "coefficient_symbols_on_common_support": ["m2_psi8", common_parameter],
            "exact_support_reduction_statement": "on this route-scoped source-action lane, the source action-role assignment witness is reduced to one still-missing coefficient-identification witness between m2_psi8 and mu_m2_plus3_segment_psi8_psi11 on the already fixed common monomial support psi8**2/2",
            "reduced_missing_object": "explicit_source_action_monomial_coefficient_identification_witness_for_m2_psi8_and_mu_m2_plus3_segment_psi8_psi11_on_common_psi8_squared_over_2_support",
            "scope": "single_direct_m2_source_action_role_common_monomial_support_only",
            "non_cancellation_clause": "R74 does not claim any global cancellation, nonzero-factor argument, or actual coefficient identification; it only exports the already shared fixed source-action monomial support",
        },
        "result": {
            "explicit_source_action_common_monomial_support_packet_present": True,
            "explicit_source_action_monomial_coefficient_identification_witness_present": False,
            "explicit_source_action_role_assignment_witness_present": False,
            "global_reduction_of_main_route_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R74",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one direct non-light m2 source action-role on the route-scoped common-monomial-support lane only",
            "boundary": "R74 does not alter the light-facing kernel channel; it only sharpens the missing source action-role assignment witness for m2_psi8 into one coefficient-identification witness on the already fixed common source-action monomial support psi8**2/2",
        },
        "classification": "explicit_source_action_common_monomial_support_packet_present_but_no_source_action_monomial_coefficient_identification_witness",
        "frontier": "R74_B1",
        "frontier_text": "on the route-scoped direct m2 sufficient lane only, the single source action-role assignment witness for m2_psi8 is sharpened into one still-missing coefficient-identification witness on the fixed common source-action monomial support psi8**2/2, while the other direct family blockers, the main route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_m2_psi8_equals_m2_psi11",
            "no_claim_that_any_common_plus3_carrier_orbit_parameter_actually_exists",
            "no_claim_that_any_source_action_coefficient_identification_witness_is_present",
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
        "stage": "R74",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI8_SOURCE_ACTION_COMMON_MONOMIAL_SUPPORT_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R74_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

