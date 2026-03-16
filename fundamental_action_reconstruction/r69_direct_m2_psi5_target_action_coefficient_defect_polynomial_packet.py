#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r69_direct_m2_psi5_target_action_coefficient_defect_polynomial_packet.json"
OUT_SUMMARY = GENERATED / "r69_direct_m2_psi5_target_action_coefficient_defect_polynomial_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r68 = load_json(
        "fundamental_action_reconstruction/generated/r68_direct_m2_psi5_target_action_common_monomial_support_packet.json"
    )

    packet = r68["direct_m2_psi5_target_action_common_monomial_support_packet"]
    target_coefficient = packet["coefficient_symbols_on_common_support"][0]
    lifted_coefficient = packet["coefficient_symbols_on_common_support"][1]
    common_support = packet["common_fixed_action_monomial_support"]
    coefficient_defect = f"({target_coefficient}) - ({lifted_coefficient})"
    action_defect = f"({coefficient_defect})*({common_support})"

    checks = [
        {
            "id": "r68_common_support_packet_present",
            "actual": r68["result"]["explicit_target_action_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R68 already exports the exact common-support packet for the target action lane",
        },
        {
            "id": "target_action_term_is_fixed",
            "actual": packet["target_action_term"],
            "expected": "m2_psi5*psi5**2/2",
            "meaning": "R68 exports the exact target action term for the attacked m2_psi5 lane",
        },
        {
            "id": "lifted_action_term_is_fixed",
            "actual": packet["lifted_action_term"],
            "expected": "mu_m2_plus3_segment_psi2_psi5*psi5**2/2",
            "meaning": "R68 exports the exact lifted action term on the same monomial support",
        },
        {
            "id": "common_support_is_fixed",
            "actual": common_support,
            "expected": "psi5**2/2",
            "meaning": "the attacked target action lane is already fixed on common support psi5**2/2",
        },
        {
            "id": "light_boundary_unchanged",
            "actual": r68["light_boundary"]["status"],
            "expected": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R68",
            "meaning": "the already closed light-facing channel remains separate from this pre-observer target action coefficient packet",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R69",
        "packet_goal": "materialize_the_exact_target_action_coefficient_defect_polynomial_whose_zero_would_give_the_missing_m2_psi5_coefficient_identification_witness_on_common_psi5_squared_over_2_support_without_dividing_by_that_support",
        "source_reports": ["R68"],
        "direct_m2_psi5_target_action_coefficient_defect_packet": {
            "slot_under_attack": "m2_psi5",
            "role_under_attack": "target_action_role",
            "target_action_term": packet["target_action_term"],
            "lifted_action_term": packet["lifted_action_term"],
            "common_fixed_action_monomial_support": common_support,
            "target_action_coefficient_symbol": target_coefficient,
            "lifted_action_coefficient_symbol": lifted_coefficient,
            "exact_coefficient_defect_polynomial": coefficient_defect,
            "exact_action_term_defect_expression": action_defect,
            "equivalent_zero_equation_on_coefficient_level": f"{coefficient_defect} = 0",
            "reduction_statement": "on this route-scoped target-action lane only, the missing coefficient-identification witness between m2_psi5 and mu_m2_plus3_segment_psi2_psi5 is reduced to one still-missing zero witness for the exact coefficient defect polynomial (m2_psi5) - (mu_m2_plus3_segment_psi2_psi5), while the common monomial support psi5**2/2 remains fixed and is not divided out",
            "reduced_missing_object": "explicit_zero_witness_for_the_direct_m2_psi5_target_action_coefficient_defect_polynomial_on_common_psi5_squared_over_2_support",
            "scope": "single_direct_m2_target_action_coefficient_defect_only",
            "non_factor_clause": "R69 does not divide by psi5**2/2, does not claim that support is nonzero, and does not claim that the coefficient defect polynomial vanishes",
        },
        "result": {
            "explicit_target_action_coefficient_defect_polynomial_present": True,
            "explicit_zero_witness_for_target_action_coefficient_defect_polynomial_present": False,
            "explicit_target_action_monomial_coefficient_identification_witness_present": False,
            "global_reduction_of_main_route_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R69",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one pre-observer non-light direct m2 target action coefficient on the already fixed common monomial support",
            "boundary": "R69 does not alter the already matched light-facing kernel channel and does not touch observer-side packets; it only materializes the exact target-action coefficient defect whose vanishing would give the missing witness for m2_psi5 on psi5**2/2",
        },
        "classification": "explicit_target_action_coefficient_defect_polynomial_present_but_no_zero_witness",
        "frontier": "R69_B1",
        "frontier_text": "on the route-scoped direct m2 sufficient lane only, the single target action-side coefficient-identification witness for m2_psi5 is sharpened into one still-missing zero witness for the exact coefficient defect polynomial on common psi5**2/2 support, while the other direct family blockers, the main route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_direct_m2_psi5_target_action_coefficient_defect_polynomial_vanishes",
            "no_claim_that_m2_psi5_equals_mu_m2_plus3_segment_psi2_psi5",
            "no_claim_that_m2_psi2_equals_m2_psi5",
            "no_claim_that_any_nonzero_factor_argument_on_psi5_squared_over_2_holds",
            "no_claim_that_the_target_eom_role_assignment_witness_is_present",
            "no_claim_that_any_other_direct_m2_pairwise_equality_holds",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R69",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI5_TARGET_ACTION_COEFFICIENT_DEFECT_POLYNOMIAL_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R69_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

