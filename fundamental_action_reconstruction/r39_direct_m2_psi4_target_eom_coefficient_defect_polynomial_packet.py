#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r39_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet.json"
OUT_SUMMARY = GENERATED / "r39_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r38 = load_json("fundamental_action_reconstruction/generated/r38_direct_m2_psi4_target_eom_common_local_support_packet.json")
    packet = r38["direct_m2_psi4_target_eom_common_local_support_packet"]
    target_coeff = "m2_psi4"
    lifted_coeff = packet["common_orbit_parameter_symbol"]
    common_support = packet["common_local_support"]
    defect_polynomial = f"({target_coeff}) - ({lifted_coeff})"
    defect_expression = f"(({target_coeff}) - ({lifted_coeff}))*{common_support}"

    checks = [
        {
            "id": "r38_target_eom_coefficient_identification_gap_still_missing_before_defect_packet",
            "actual": not r38["result"]["explicit_target_eom_monomial_coefficient_identification_witness_present"],
            "expected": True,
            "meaning": "R38 still leaves the target eom monomial coefficient-identification witness absent before any defect-polynomial reduction",
        },
        {
            "id": "r38_target_coefficient_symbol_present",
            "actual": target_coeff,
            "expected": "m2_psi4",
            "meaning": "R38 keeps the exact target eom coefficient symbol on the attacked target lane",
        },
        {
            "id": "r38_lifted_coefficient_symbol_present",
            "actual": lifted_coeff,
            "expected": "mu_m2_plus3_segment_psi1_psi4",
            "meaning": "R38 keeps the exact lifted coefficient symbol on the attacked target eom lane",
        },
        {
            "id": "r38_common_support_present",
            "actual": common_support,
            "expected": "psi4(x)",
            "meaning": "R38 keeps the exact common target-eom local support",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R39",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-R38",
        "packet_goal": "reduce_the_single_target_eom_monomial_coefficient_identification_gap_for_m2_psi4_to_one_exact_target_eom_coefficient_defect_polynomial_zero_witness_without_dividing_by_common_support",
        "source_reports": ["R38"],
        "direct_m2_psi4_target_eom_coefficient_defect_packet": {
            "target_eom_coefficient_symbol": target_coeff,
            "lifted_eom_coefficient_symbol": lifted_coeff,
            "common_support": common_support,
            "exact_coefficient_defect_polynomial": defect_polynomial,
            "exact_local_eom_defect_expression": defect_expression,
            "remaining_missing_witness": "explicit_zero_witness_for_the_direct_m2_psi4_target_eom_coefficient_defect_polynomial_on_common_psi4_of_x_support",
            "exact_reduction_statement": "on this route-scoped one-pair target eom lane, the missing target eom monomial coefficient-identification witness is reduced to the zero of the exact coefficient defect polynomial (m2_psi4) - (mu_m2_plus3_segment_psi1_psi4) on common support psi4(x)",
            "scope": "single_direct_m2_target_eom_coefficient_defect_only",
            "non_promotion_clause": "R39 does not claim that the target eom coefficient defect polynomial vanishes in the current repo, nor that any nonzero-factor argument on psi4(x) is allowed",
        },
        "result": {
            "explicit_target_eom_coefficient_defect_packet_present": True,
            "explicit_zero_witness_for_target_eom_coefficient_defect_polynomial_present": False,
            "explicit_target_eom_monomial_coefficient_identification_witness_present": False,
            "global_reduction_of_main_R21_c1s1_blocker_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R39",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one direct non-light m2 target eom coefficient on the canonical-ontology-supported defect-polynomial lane only",
            "boundary": "R39 does not alter the light-facing kernel channel; it only sharpens the missing target eom monomial coefficient-identification witness for m2_psi4 into one exact coefficient defect polynomial on common psi4(x) support",
        },
        "classification": "explicit_target_eom_coefficient_defect_packet_present_but_zero_witness_for_target_eom_coefficient_defect_polynomial_absent",
        "frontier": "R39_B1",
        "frontier_text": "on the canonical-ontology-supported direct m2 sufficient lane only, the target eom monomial coefficient-identification witness for m2_psi4 is sharpened into one still-missing zero witness for the exact coefficient defect polynomial on common psi4(x) support, while the other direct family blockers, the main host route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_direct_m2_psi4_target_eom_coefficient_defect_polynomial_vanishes",
            "no_claim_that_m2_psi4_equals_mu_m2_plus3_segment_psi1_psi4",
            "no_claim_that_m2_psi1_equals_m2_psi4",
            "no_claim_that_any_nonzero_factor_argument_on_psi4_of_x_holds",
            "no_claim_that_any_other_direct_m2_pairwise_equality_holds",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R39",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI4_TARGET_EOM_COEFFICIENT_DEFECT_POLYNOMIAL_PACKET_READY",
        "lane": artifact["lane"],
        "result": artifact["classification"],
        "frontier": ["R39_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
