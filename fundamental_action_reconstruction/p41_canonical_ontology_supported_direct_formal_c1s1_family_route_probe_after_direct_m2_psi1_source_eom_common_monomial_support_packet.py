#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p41_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_source_eom_common_monomial_support_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p41_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_source_eom_common_monomial_support_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p40 = load_json(
        "fundamental_action_reconstruction/generated/p40_axiom_lane_direct_formal_c1s1_family_route_probe_after_preobserver_source_action_coherence_instance.json"
    )
    r33 = load_json("fundamental_action_reconstruction/generated/r33_direct_m2_psi1_source_eom_common_monomial_support_packet.json")

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_source_eom_monomial_coefficient_identification_witness_for_m2_psi1_and_mu_m2_plus3_segment_psi1_psi4_on_common_psi1_of_x_support",
        "explicit_assignment_witness_of_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
        "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10",
        "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
        "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    checks = [
        {
            "id": "p40_local_source_action_closure_present",
            "actual": p40["route_state"]["attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only"],
            "expected": True,
            "meaning": "P40 already records the local closure of the attacked source-action blocker on the canonical-ontology-supported lane only",
        },
        {
            "id": "p40_route_still_not_closed",
            "actual": p40["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_ONLY_ATTACKED_SOURCE_ACTION_BLOCKER_CLOSED_ROUTE_STILL_NOT_CLOSED_AFTER_AX10",
            "meaning": "before R33, the canonical-ontology-supported route still remained non-closed",
        },
        {
            "id": "r33_source_eom_common_support_packet_present",
            "actual": r33["result"]["explicit_source_eom_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R33 adds the exact common-support packet for the attacked source eom-role witness",
        },
        {
            "id": "r33_source_eom_coefficient_identification_witness_still_absent",
            "actual": r33["result"]["explicit_source_eom_monomial_coefficient_identification_witness_present"],
            "expected": False,
            "meaning": "R33 does not claim that the source eom-side coefficient-identification witness already holds",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "P41",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-and-R33",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_reducing_the_attacked_source_eom_role_assignment_witness_to_one_common_local_eom_support_coefficient_gap",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_ONLY_ATTACKED_SOURCE_ACTION_BLOCKER_CLOSED_AND_SOURCE_EOM_GAP_REDUCED_TO_COMMON_PSI1_OF_X_SUPPORT_ROUTE_STILL_NOT_CLOSED_AFTER_R33",
        "reason": "AX10 closes the attacked source-action blocker on the canonical-ontology-supported pre-observer lane only, and R33 reduces the attacked source eom-role witness to one coefficient-identification witness on common psi1(x) support, but the source eom-side coefficient-identification witness remains absent, the target-side assignment witness remains absent, the direct g4, g6, and gY family zero witnesses remain absent, the other direct m2 pairwise witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "source_eom_gap_reduced_to_common_local_support": True,
            "source_eom_monomial_coefficient_identification_witness_present": False,
            "target_slot_assignment_witness_present": False,
            "direct_g4_family_zero_witness_present": False,
            "direct_g6_family_zero_witness_present": False,
            "direct_gY_family_zero_witness_present": False,
            "direct_m2_psi7_psi10_pairwise_witness_present": False,
            "direct_m2_psi2_psi5_pairwise_witness_present": False,
            "direct_m2_psi8_psi11_pairwise_witness_present": False,
            "pair1_c1c1_zero_witness_present": False,
            "pair1_s1s1_zero_witness_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "decomposition_of_P40_source_eom_gap": {
            "from_P40": "explicit_source_eom_role_assignment_witness_for_m2_psi1_to_mu_m2_plus3_segment_psi1_psi4",
            "reduced_by_R33_into": "explicit_source_eom_monomial_coefficient_identification_witness_for_m2_psi1_and_mu_m2_plus3_segment_psi1_psi4_on_common_psi1_of_x_support",
            "scope": "single_direct_m2_source_eom_role_common_local_support_only",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P41",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "decomposition_of_P40_source_eom_gap": artifact["decomposition_of_P40_source_eom_gap"],
        "remaining_missing_upstream_objects": remaining_missing,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
