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
    / "p45_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_action_common_monomial_support_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p45_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_action_common_monomial_support_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p44 = load_json(
        "fundamental_action_reconstruction/generated/p44_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_common_plus3_assignment_role_split_packet.json"
    )
    r36 = load_json("fundamental_action_reconstruction/generated/r36_direct_m2_psi4_target_action_common_monomial_support_packet.json")

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_target_action_monomial_coefficient_identification_witness_for_m2_psi4_and_mu_m2_plus3_segment_psi1_psi4_on_common_psi4_squared_over_2_support",
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
            "id": "p44_target_action_role_gap_was_still_missing",
            "actual": "explicit_target_action_role_assignment_witness_for_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4"
            in p44["remaining_missing_objects"],
            "expected": True,
            "meaning": "P44 still had the target action-role assignment witness for m2_psi4 as a missing object",
        },
        {
            "id": "r36_target_action_common_support_packet_present",
            "actual": r36["result"]["explicit_target_action_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R36 adds the exact common target-action support packet for the attacked target-side lane",
        },
        {
            "id": "r36_target_action_monomial_coefficient_identification_witness_still_absent",
            "actual": r36["result"]["explicit_target_action_monomial_coefficient_identification_witness_present"],
            "expected": False,
            "meaning": "R36 does not claim that the target action-side coefficient-identification witness holds",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "P45",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-AX11-R35-and-R36",
        "goal": "rerun_the_canonical_ontology_supported_direct_family_route_after_exact_target_action_common_support_reduction_of_the_attacked_target_action_role_assignment_witness_for_m2_psi4",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_CLOSED_AND_TARGET_ACTION_GAP_REDUCED_TO_COMMON_PSI4_SQUARED_OVER_2_SUPPORT_ROUTE_STILL_NOT_CLOSED_AFTER_R36",
        "reason": "AX10 and AX11 locally close the attacked source-side blockers on the canonical-ontology-supported pre-observer lane, R35 splits the target-slot assignment witness for m2_psi4 into target action-role and target eom-role assignment witnesses, and R36 reduces only the target action-role witness to one target action monomial coefficient-identification gap on common psi4**2/2 support; the target eom-role assignment witness remains absent, the direct g4, g6, and gY family zero witnesses remain absent, the other direct m2 pairwise witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_source_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "target_action_monomial_coefficient_identification_witness_present": False,
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
        "decomposition_of_P44_target_action_gap": {
            "from_P44": "explicit_target_action_role_assignment_witness_for_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
            "reduced_by_R36_into": "explicit_target_action_monomial_coefficient_identification_witness_for_m2_psi4_and_mu_m2_plus3_segment_psi1_psi4_on_common_psi4_squared_over_2_support",
            "scope": "single_direct_m2_target_action_role_common_support_only",
        },
        "remaining_missing_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P45",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "decomposition_of_P44_target_action_gap": artifact["decomposition_of_P44_target_action_gap"],
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
