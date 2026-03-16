#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r77_direct_m2_psi8_source_eom_coefficient_defect_polynomial_packet.json"
OUT_SUMMARY = GENERATED / "r77_direct_m2_psi8_source_eom_coefficient_defect_polynomial_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r76 = load_json("fundamental_action_reconstruction/generated/r76_direct_m2_psi8_source_eom_common_monomial_support_packet.json")

    packet = r76["direct_m2_psi8_source_eom_common_monomial_support_packet"]
    source_term = packet["source_eom_term"]
    lifted_term = packet["declared_lifted_eom_term"]
    source_symbol = "m2_psi8"
    lifted_symbol = packet["common_parameter_symbol"]
    common_support = packet["fixed_common_local_eom_support"]
    defect_polynomial = f"({source_symbol}) - ({lifted_symbol})"
    defect_expression = f"({defect_polynomial})*{common_support}"

    checks = [
        {
            "id": "r76_source_term_present",
            "actual": source_term,
            "expected": "m2_psi8*psi8(x)",
            "meaning": "R76 exports the exact source local eom term on the attacked route",
        },
        {
            "id": "r76_lifted_term_present",
            "actual": lifted_term,
            "expected": "mu_m2_plus3_segment_psi8_psi11*psi8(x)",
            "meaning": "R76 exports the exact lifted local eom term on the same route",
        },
        {
            "id": "r76_common_support_present",
            "actual": common_support,
            "expected": "psi8(x)",
            "meaning": "R76 exports the exact common local support psi8(x)",
        },
        {
            "id": "r76_source_eom_coefficient_identification_witness_still_absent",
            "actual": r76["result"]["explicit_source_eom_monomial_coefficient_identification_witness_present"],
            "expected": False,
            "meaning": "R76 does not claim that the source eom-side coefficient-identification witness already holds",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R77",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-R72",
        "packet_goal": "reduce_the_single_source_eom_coefficient_identification_gap_to_one_exact_coefficient_defect_polynomial_without_dividing_by_the_local_support_and_without_any_nonzero_factor_claim",
        "source_reports": ["R76"],
        "direct_m2_psi8_source_eom_coefficient_defect_packet": {
            "source_eom_coefficient_symbol": source_symbol,
            "lifted_eom_coefficient_symbol": lifted_symbol,
            "common_fixed_local_support": common_support,
            "exact_source_eom_term": source_term,
            "exact_lifted_eom_term": lifted_term,
            "exact_coefficient_defect_polynomial": defect_polynomial,
            "exact_local_eom_defect_expression": defect_expression,
            "exact_zero_equation_needed_for_missing_witness": f"{defect_polynomial} = 0",
            "scope": "single_direct_m2_source_eom_coefficient_gap_only",
            "non_promotion_clause": "R77 does not claim that the source-eom defect polynomial vanishes, nor that psi8(x) may be divided out or treated as a nonzero factor",
        },
        "result": {
            "explicit_source_eom_coefficient_defect_polynomial_present": True,
            "explicit_zero_witness_for_source_eom_coefficient_defect_polynomial_present": False,
            "explicit_source_eom_monomial_coefficient_identification_witness_present": False,
            "strict_core_promotion": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R77",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one pre-observer non-light source eom coefficient on the route-scoped source-eom lane only",
            "boundary": "R77 keeps light before observer explicit and sharpens only the source eom-side coefficient gap without any factor or division claim",
        },
        "classification": "explicit_source_eom_coefficient_defect_polynomial_present_but_zero_witness_absent",
        "frontier": "R77_B1",
        "frontier_text": "on the canonical-ontology-supported direct route after R72, the single source eom-side coefficient-identification gap for m2_psi8 is sharpened into one still-missing zero witness for an exact defect polynomial on common psi8(x) support, while the target-side, other direct-family blockers, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_source_eom_defect_polynomial_vanishes",
            "no_claim_that_m2_psi8_equals_mu_m2_plus3_segment_psi8_psi11",
            "no_claim_that_m2_psi8_equals_m2_psi11",
            "no_claim_that_any_field_division_or_nonzero_factor_argument_on_psi8_of_x_holds",
            "no_claim_that_the_target_side_assignment_witness_is_present",
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
        "stage": "R77",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI8_SOURCE_EOM_COEFFICIENT_DEFECT_POLYNOMIAL_PACKET_READY",
        "lane": artifact["lane"],
        "result": artifact["classification"],
        "frontier": ["R77_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

