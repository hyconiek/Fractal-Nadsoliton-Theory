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
    / "p44_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_common_plus3_assignment_role_split_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p44_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_common_plus3_assignment_role_split_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p43 = load_json(
        "fundamental_action_reconstruction/generated/p43_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_source_eom_coherence_instance.json"
    )
    r35 = load_json("fundamental_action_reconstruction/generated/r35_direct_m2_psi4_common_plus3_assignment_role_split_packet.json")

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_target_action_role_assignment_witness_for_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
        "explicit_target_eom_role_assignment_witness_for_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
        "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10",
        "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
        "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    checks = [
        {
            "id": "p43_target_slot_assignment_gap_was_still_missing",
            "actual": "explicit_assignment_witness_of_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4" in p43["remaining_missing_objects"],
            "expected": True,
            "meaning": "P43 still had the target-slot assignment witness for m2_psi4 as a missing object",
        },
        {
            "id": "r35_target_assignment_role_split_packet_present",
            "actual": r35["result"]["explicit_target_slot_assignment_role_split_packet_present"],
            "expected": True,
            "meaning": "R35 adds the exact role split of the target-slot assignment witness for m2_psi4",
        },
        {
            "id": "r35_target_action_role_assignment_witness_still_absent",
            "actual": r35["result"]["explicit_target_action_role_assignment_witness_present"],
            "expected": False,
            "meaning": "R35 does not claim that the target action-role assignment witness holds",
        },
        {
            "id": "r35_target_eom_role_assignment_witness_still_absent",
            "actual": r35["result"]["explicit_target_eom_role_assignment_witness_present"],
            "expected": False,
            "meaning": "R35 does not claim that the target eom-role assignment witness holds",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "P44",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-AX11-and-R35",
        "goal": "rerun_the_canonical_ontology_supported_direct_family_route_after_exact_target_role_splitting_of_the_target_slot_assignment_witness_for_m2_psi4",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_CLOSED_BUT_TARGET_M2_PSI4_ROLE_ASSIGNMENT_GAPS_REMAIN_AFTER_R35",
        "reason": "AX10 and AX11 locally close the attacked source-side blockers on the canonical-ontology-supported pre-observer lane, but R35 only splits the target-slot assignment witness for m2_psi4 into target action-role and target eom-role assignment witnesses, and both remain absent; the direct g4, g6, and gY family zero witnesses remain absent, the other direct m2 pairwise witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_source_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "target_action_role_assignment_witness_present": False,
            "target_eom_role_assignment_witness_present": False,
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
        "decomposition_of_P43_target_gap": {
            "from_P43": "explicit_assignment_witness_of_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
            "reduced_by_R35_into": [
                "explicit_target_action_role_assignment_witness_for_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
                "explicit_target_eom_role_assignment_witness_for_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
            ],
            "scope": "single_direct_m2_target_slot_m2_psi4_common_plus3_assignment_role_split_only",
        },
        "remaining_missing_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P44",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "decomposition_of_P43_target_gap": artifact["decomposition_of_P43_target_gap"],
        "remaining_missing_objects": remaining_missing,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
