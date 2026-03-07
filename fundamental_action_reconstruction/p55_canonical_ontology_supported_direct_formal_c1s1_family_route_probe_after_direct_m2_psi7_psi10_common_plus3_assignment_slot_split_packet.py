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
    / "p55_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi7_psi10_common_plus3_assignment_slot_split_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p55_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi7_psi10_common_plus3_assignment_slot_split_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p54 = load_json(
        "fundamental_action_reconstruction/generated/p54_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi7_psi10_common_plus3_orbit_parameter_sufficient_route_packet.json"
    )
    r43 = load_json(
        "fundamental_action_reconstruction/generated/r43_direct_m2_psi7_psi10_common_plus3_assignment_slot_split_packet.json"
    )

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_assignment_witness_of_m2_psi7_to_mu_m2_plus3_segment_psi7_psi10",
        "explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10",
        "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
        "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    checks = [
        {
            "id": "p54_combined_assignment_gap_was_still_missing",
            "actual": p54["decomposition_of_remaining_pairwise_gap"]["reduced_by_R42_into"],
            "expected": "explicit_assignment_witness_of_m2_psi7_and_m2_psi10_to_one_common_plus3_carrier_segment_parameter",
            "meaning": "before R43, P54 still tracked the combined assignment witness for m2_psi7 / m2_psi10 as the narrowest local gap",
        },
        {
            "id": "r43_assignment_slot_split_packet_present",
            "actual": r43["result"]["explicit_common_plus3_assignment_slot_split_packet_present"],
            "expected": True,
            "meaning": "R43 adds the exact slotwise split of the one-pair combined assignment witness",
        },
        {
            "id": "r43_source_slot_assignment_witness_still_absent",
            "actual": r43["result"]["explicit_assignment_witness_of_m2_psi7_to_mu_m2_plus3_segment_psi7_psi10_present"],
            "expected": False,
            "meaning": "R43 does not claim that the source-slot assignment witness holds",
        },
        {
            "id": "r43_target_slot_assignment_witness_still_absent",
            "actual": r43["result"]["explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10_present"],
            "expected": False,
            "meaning": "R43 does not claim that the target-slot assignment witness holds",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "P55",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-AX11-AX12-AX13-R40-R41-R42-and-R43",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_R43_direct_m2_psi7_psi10_common_plus3_assignment_slot_split_packet",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_DIRECT_M2_PSI7_PSI10_PAIRWISE_GAP_REDUCED_TO_SLOTWISE_ASSIGNMENT_PACKET_ROUTE_STILL_NOT_CLOSED_AFTER_R43",
        "reason": "AX10, AX11, AX12, and AX13 close the attacked source-side and attacked m2_psi4 target-side blockers on the canonical-ontology-supported pre-observer lane only, R40 reduces the direct pairwise target m2_psi7 = m2_psi10 to one coefficient-identification witness, R41 reduces that witness to one common-source or symbol-identification witness, R42 reduces that witness to one combined assignment witness to a common plus3 carrier-segment parameter, and R43 reduces that combined witness further to two slotwise assignment witnesses, but both slotwise witnesses remain absent, the two other direct m2 pairwise witnesses remain absent, the direct g4, g6, and gY family zero witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_source_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "direct_m2_psi7_psi10_assignment_slot_split_packet_present": True,
            "direct_m2_psi7_assignment_witness_present": False,
            "direct_m2_psi10_assignment_witness_present": False,
            "direct_m2_psi2_psi5_pairwise_witness_present": False,
            "direct_m2_psi8_psi11_pairwise_witness_present": False,
            "direct_g4_family_zero_witness_present": False,
            "direct_g6_family_zero_witness_present": False,
            "direct_gY_family_zero_witness_present": False,
            "pair1_c1c1_zero_witness_present": False,
            "pair1_s1s1_zero_witness_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "decomposition_of_remaining_pairwise_gap": {
            "from_P54": "explicit_assignment_witness_of_m2_psi7_and_m2_psi10_to_one_common_plus3_carrier_segment_parameter",
            "reduced_by_R43_into": [
                "explicit_assignment_witness_of_m2_psi7_to_mu_m2_plus3_segment_psi7_psi10",
                "explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10",
            ],
            "scope": "single_direct_m2_pairwise_gap_only",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P55",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "decomposition_of_remaining_pairwise_gap": artifact["decomposition_of_remaining_pairwise_gap"],
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
