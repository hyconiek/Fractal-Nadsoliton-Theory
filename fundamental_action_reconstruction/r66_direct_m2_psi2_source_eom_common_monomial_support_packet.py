#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r66_direct_m2_psi2_source_eom_common_monomial_support_packet.json"
OUT_SUMMARY = GENERATED / "r66_direct_m2_psi2_source_eom_common_monomial_support_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r62 = load_json("fundamental_action_reconstruction/generated/r62_direct_m2_psi2_common_plus3_assignment_role_split_packet.json")

    role_packet = r62["direct_m2_psi2_common_plus3_assignment_role_split_packet"]
    source_eom_term = role_packet["source_role_support"]["canonical_eom_source_term"]
    lifted_eom_term = role_packet["declared_rolewise_parameter_lifts"]["eom_role_lifted_term"]
    common_parameter = role_packet["common_orbit_parameter_symbol"]
    common_support = "psi2(x)"
    missing_witness = (
        "explicit_source_eom_monomial_coefficient_identification_witness_for_"
        "m2_psi2_and_mu_m2_plus3_segment_psi2_psi5_on_common_psi2_of_x_support"
    )

    checks = [
        {
            "id": "r62_source_eom_role_term_present",
            "actual": source_eom_term,
            "expected": "m2_psi2*psi2(x)",
            "meaning": "R62 exports the exact source local eom term for the attacked m2_psi2 role",
        },
        {
            "id": "r62_lifted_source_eom_term_present",
            "actual": lifted_eom_term,
            "expected": "mu_m2_plus3_segment_psi2_psi5*psi2(x)",
            "meaning": "R62 exports the exact declared lifted local eom term on the same route",
        },
        {
            "id": "r62_source_eom_role_assignment_witness_still_absent",
            "actual": r62["result"]["explicit_source_eom_role_assignment_witness_present"],
            "expected": False,
            "meaning": "R62 does not claim that the source eom-role assignment witness already holds",
        },
        {
            "id": "r62_common_parameter_symbol_present",
            "actual": common_parameter,
            "expected": "mu_m2_plus3_segment_psi2_psi5",
            "meaning": "R62 preserves the declared common plus3 carrier-segment parameter on the one-pair route",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R66",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-R62",
        "packet_goal": "materialize_the_exact_common_local_eom_support_packet_for_the_single_attacked_source_eom_role_assignment_witness_without_claiming_the_assignment_or_any_field_division_argument",
        "source_reports": ["R62"],
        "direct_m2_psi2_source_eom_common_monomial_support_packet": {
            "source_eom_role_under_attack": "explicit_source_eom_role_assignment_witness_for_m2_psi2_to_mu_m2_plus3_segment_psi2_psi5",
            "source_eom_term": source_eom_term,
            "declared_lifted_eom_term": lifted_eom_term,
            "common_parameter_symbol": common_parameter,
            "fixed_common_local_eom_support": common_support,
            "reduced_missing_witness": missing_witness,
            "scope": "single_direct_m2_source_eom_role_common_local_support_only",
            "non_promotion_clause": "R66 does not claim that the source eom-role assignment witness holds, nor that the fixed local eom support may be divided out or treated as nonzero",
        },
        "result": {
            "explicit_source_eom_common_monomial_support_packet_present": True,
            "explicit_source_eom_monomial_coefficient_identification_witness_present": False,
            "explicit_source_eom_role_assignment_witness_present": False,
            "strict_core_promotion": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R66",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one pre-observer non-light direct m2 source eom role on the route-scoped source-eom common-support lane only",
            "boundary": "R66 keeps light before observer explicit and sharpens only the source eom-side gap without any factor or division claim",
        },
        "classification": "explicit_source_eom_common_monomial_support_packet_present_but_source_eom_monomial_coefficient_identification_witness_absent",
        "frontier": "R66_B1",
        "frontier_text": "on the canonical-ontology-supported direct route after R62, the single source eom-role witness for m2_psi2 is sharpened into one still-missing coefficient-identification witness on common psi2(x) support, while the target-side, other direct-family blockers, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_m2_psi2_equals_m2_psi5",
            "no_claim_that_any_common_plus3_carrier_segment_parameter_actually_exists",
            "no_claim_that_the_source_eom_side_coefficient_identification_witness_is_present",
            "no_claim_that_any_field_division_or_nonzero_factor_argument_on_psi2_of_x_holds",
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
        "stage": "R66",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI2_SOURCE_EOM_COMMON_MONOMIAL_SUPPORT_PACKET_READY",
        "lane": artifact["lane"],
        "result": artifact["classification"],
        "frontier": ["R66_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

