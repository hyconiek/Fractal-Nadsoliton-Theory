#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r71_direct_m2_psi5_target_eom_coefficient_defect_polynomial_packet.json"
OUT_SUMMARY = GENERATED / "r71_direct_m2_psi5_target_eom_coefficient_defect_polynomial_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r70 = load_json("fundamental_action_reconstruction/generated/r70_direct_m2_psi5_target_eom_common_monomial_support_packet.json")

    packet = r70["direct_m2_psi5_target_eom_common_monomial_support_packet"]
    target_eom_term = packet["target_eom_term"]
    lifted_eom_term = packet["declared_lifted_eom_term"]
    common_support = packet["fixed_common_local_eom_support"]

    target_coefficient = "m2_psi5"
    lifted_coefficient = packet["common_parameter_symbol"]
    coefficient_defect = f"({target_coefficient}) - ({lifted_coefficient})"
    eom_defect = f"({coefficient_defect})*{common_support}"

    checks = [
        {
            "id": "r70_common_support_packet_present",
            "actual": r70["result"]["explicit_target_eom_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R70 already exports the exact common-support packet for the target eom lane",
        },
        {
            "id": "target_eom_term_is_fixed",
            "actual": target_eom_term,
            "expected": "m2_psi5*psi5(x)",
            "meaning": "R70 exports the exact target local eom term for the attacked m2_psi5 lane",
        },
        {
            "id": "lifted_eom_term_is_fixed",
            "actual": lifted_eom_term,
            "expected": "mu_m2_plus3_segment_psi2_psi5*psi5(x)",
            "meaning": "R70 exports the exact lifted local eom term on the same support",
        },
        {
            "id": "common_support_is_fixed",
            "actual": common_support,
            "expected": "psi5(x)",
            "meaning": "the attacked target eom lane is already fixed on common support psi5(x)",
        },
        {
            "id": "light_boundary_unchanged",
            "actual": r70["light_boundary"]["status"],
            "expected": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R70",
            "meaning": "the already closed light-facing channel remains separate from this pre-observer target eom coefficient packet",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R71",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-R63",
        "packet_goal": "reduce_the_single_target_eom_coefficient_identification_gap_to_one_exact_coefficient_defect_polynomial_without_dividing_by_the_local_support_and_without_any_nonzero_factor_claim",
        "source_reports": ["R70"],
        "direct_m2_psi5_target_eom_coefficient_defect_packet": {
            "target_eom_coefficient_symbol": target_coefficient,
            "lifted_eom_coefficient_symbol": lifted_coefficient,
            "common_fixed_local_support": common_support,
            "exact_target_eom_term": target_eom_term,
            "exact_lifted_eom_term": lifted_eom_term,
            "exact_coefficient_defect_polynomial": coefficient_defect,
            "exact_local_eom_defect_expression": eom_defect,
            "exact_zero_equation_needed_for_missing_witness": f"{coefficient_defect} = 0",
            "scope": "single_direct_m2_target_eom_coefficient_defect_only",
            "non_promotion_clause": "R71 does not claim that the target-eom defect polynomial vanishes, nor that psi5(x) may be divided out or treated as a nonzero factor",
        },
        "result": {
            "explicit_target_eom_coefficient_defect_polynomial_present": True,
            "explicit_zero_witness_for_target_eom_coefficient_defect_polynomial_present": False,
            "explicit_target_eom_monomial_coefficient_identification_witness_present": False,
            "strict_core_promotion": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R71",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one pre-observer non-light direct m2 target eom coefficient on the already fixed common local support",
            "boundary": "R71 keeps light before observer explicit and sharpens only the target eom-side coefficient gap without any factor or division claim",
        },
        "classification": "explicit_target_eom_coefficient_defect_polynomial_present_but_zero_witness_absent",
        "frontier": "R71_B1",
        "frontier_text": "on the canonical-ontology-supported direct route after R63, the single target eom-side coefficient-identification gap for m2_psi5 is sharpened into one still-missing zero witness for an exact defect polynomial on common psi5(x) support, while the source-side, other direct-family blockers, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_target_eom_defect_polynomial_vanishes",
            "no_claim_that_m2_psi5_equals_mu_m2_plus3_segment_psi2_psi5",
            "no_claim_that_m2_psi2_equals_m2_psi5",
            "no_claim_that_any_field_division_or_nonzero_factor_argument_on_psi5_of_x_holds",
            "no_claim_that_any_other_direct_m2_pairwise_equality_holds",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
            "no_claim_that_ToE_is_closed",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R71",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI5_TARGET_EOM_COEFFICIENT_DEFECT_POLYNOMIAL_PACKET_READY",
        "lane": artifact["lane"],
        "result": artifact["classification"],
        "frontier": ["R71_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

