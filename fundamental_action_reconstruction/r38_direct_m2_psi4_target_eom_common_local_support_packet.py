#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r38_direct_m2_psi4_target_eom_common_local_support_packet.json"
OUT_SUMMARY = GENERATED / "r38_direct_m2_psi4_target_eom_common_local_support_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r26 = load_json("fundamental_action_reconstruction/generated/r26_direct_m2_psi1_psi4_role_matching_packet.json")
    r35 = load_json("fundamental_action_reconstruction/generated/r35_direct_m2_psi4_common_plus3_assignment_role_split_packet.json")
    ax12 = load_json("fundamental_action_reconstruction/generated/ax12_canonical_ontology_supported_preobserver_target_action_coherence_instance.json")

    role_packet = r26["direct_m2_psi1_psi4_role_matching_packet"]
    target_eom_term = role_packet["canonical_eom_role_match"]["target_term"]
    lifted_eom_term = r35["direct_m2_psi4_common_plus3_assignment_role_split_packet"][
        "declared_rolewise_parameter_lifts"
    ]["eom_role_lifted_term"]
    common_parameter = r35["direct_m2_psi4_common_plus3_assignment_role_split_packet"]["common_orbit_parameter_symbol"]

    checks = [
        {
            "id": "ax12_attacked_target_action_blocker_already_closed",
            "actual": ax12["result"]["canonical_ontology_supported_target_action_coefficient_defect_zero_witness_present"],
            "expected": True,
            "meaning": "AX12 already closes the attacked target-action blocker on the canonical-ontology-supported lane",
        },
        {
            "id": "r35_target_eom_role_assignment_witness_still_absent_before_common_support_reduction",
            "actual": r35["result"]["explicit_target_eom_role_assignment_witness_present"],
            "expected": False,
            "meaning": "R35 still leaves the target eom-role assignment witness absent before any finer common-support reduction",
        },
        {
            "id": "r26_target_eom_role_present",
            "actual": target_eom_term,
            "expected": "m2_psi4*psi4(x)",
            "meaning": "R26 exports the exact target local eom role for m2_psi4",
        },
        {
            "id": "r35_declared_target_eom_lift_present",
            "actual": lifted_eom_term,
            "expected": "mu_m2_plus3_segment_psi1_psi4*psi4(x)",
            "meaning": "R35 exports the exact declared lifted target eom term for the attacked target-side lane",
        },
        {
            "id": "r35_common_parameter_symbol_present",
            "actual": common_parameter,
            "expected": "mu_m2_plus3_segment_psi1_psi4",
            "meaning": "R35 preserves the declared common plus3 carrier-segment parameter symbol on the one-pair target lane",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R38",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX12",
        "packet_goal": "reduce_the_single_target_eom_role_assignment_witness_for_m2_psi4_to_one_target_eom_monomial_coefficient_identification_witness_on_fixed_common_psi4_of_x_support_without_claiming_that_witness",
        "source_reports": ["AX12", "R26", "R35"],
        "direct_m2_psi4_target_eom_common_local_support_packet": {
            "slot_under_attack": "m2_psi4",
            "target_eom_role_assignment_witness_under_attack": "explicit_target_eom_role_assignment_witness_for_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
            "target_eom_term": target_eom_term,
            "declared_lifted_target_eom_term": lifted_eom_term,
            "common_orbit_parameter_symbol": common_parameter,
            "common_local_support": "psi4(x)",
            "remaining_missing_witness": "explicit_target_eom_monomial_coefficient_identification_witness_for_m2_psi4_and_mu_m2_plus3_segment_psi1_psi4_on_common_psi4_of_x_support",
            "exact_reduction_statement": "on this route-scoped one-pair target lane, the target eom-role assignment witness can only be resolved by identifying the coefficient labels m2_psi4 and mu_m2_plus3_segment_psi1_psi4 on the already fixed common target local eom support psi4(x)",
            "scope": "single_direct_m2_target_eom_role_common_support_only",
            "non_promotion_clause": "R38 does not claim that the target eom-side coefficient-identification witness is present in the current repo, nor that any field-division or nonzero-factor argument on psi4(x) holds",
        },
        "result": {
            "explicit_target_eom_common_local_support_packet_present": True,
            "explicit_target_eom_monomial_coefficient_identification_witness_present": False,
            "explicit_target_eom_role_assignment_witness_present": False,
            "global_reduction_of_main_R21_c1s1_blocker_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R38",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one direct non-light m2 target eom-role on the canonical-ontology-supported common-support lane only",
            "boundary": "R38 does not alter the light-facing kernel channel; it only sharpens the missing target eom-role assignment witness for m2_psi4 into a target eom monomial coefficient-identification witness on common psi4(x) support",
        },
        "classification": "explicit_target_eom_common_local_support_packet_present_but_target_eom_monomial_coefficient_identification_witness_absent",
        "frontier": "R38_B1",
        "frontier_text": "on the canonical-ontology-supported direct m2 sufficient lane only, the target eom-role assignment witness for m2_psi4 is sharpened into one still-missing target eom monomial coefficient-identification witness on common psi4(x) support, while the other direct family blockers, the main host route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_m2_psi1_equals_m2_psi4",
            "no_claim_that_any_common_plus3_carrier_orbit_parameter_actually_exists",
            "no_claim_that_the_target_eom_side_coefficient_identification_witness_is_present",
            "no_claim_that_any_field_division_or_nonzero_factor_argument_holds",
            "no_claim_that_any_other_direct_m2_pairwise_equality_holds",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R38",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI4_TARGET_EOM_COMMON_LOCAL_SUPPORT_PACKET_READY",
        "lane": artifact["lane"],
        "result": artifact["classification"],
        "frontier": ["R38_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
