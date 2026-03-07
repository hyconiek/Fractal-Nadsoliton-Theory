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
    / "p60_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi7_source_eom_common_monomial_support_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p60_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi7_source_eom_common_monomial_support_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p59 = load_json(
        "fundamental_action_reconstruction/generated/p59_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi7_source_action_coherence_instance.json"
    )
    r47 = load_json(
        "fundamental_action_reconstruction/generated/r47_direct_m2_psi7_source_eom_common_monomial_support_packet.json"
    )

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_source_eom_monomial_coefficient_identification_witness_for_m2_psi7_and_mu_m2_plus3_segment_psi7_psi10_on_common_psi7_of_x_support",
        "explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10",
        "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
        "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    checks = [
        {
            "id": "p59_source_eom_role_assignment_gap_was_still_missing",
            "actual": "explicit_source_eom_role_assignment_witness_for_m2_psi7_to_mu_m2_plus3_segment_psi7_psi10" in p59["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P59 still had the source eom-role assignment witness for m2_psi7 as a missing object",
        },
        {
            "id": "r47_source_eom_common_support_packet_present",
            "actual": r47["result"]["explicit_source_eom_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R47 adds the exact common local support packet for the attacked source eom lane",
        },
        {
            "id": "r47_source_eom_monomial_coefficient_identification_witness_still_absent",
            "actual": r47["result"]["explicit_source_eom_monomial_coefficient_identification_witness_present"],
            "expected": False,
            "meaning": "R47 does not claim that the source eom coefficient-identification witness is already present",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "P60",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-AX11-AX12-AX13-R40-R41-R42-R43-R44-R45-R46-AX14-and-R47",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_R47_direct_m2_psi7_source_eom_common_monomial_support_packet",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_DIRECT_M2_PSI7_SOURCE_EOM_GAP_REDUCED_TO_COMMON_MONOMIAL_SUPPORT_PACKET_ROUTE_STILL_NOT_CLOSED_AFTER_R47",
        "reason": "AX10, AX11, AX12, AX13, and AX14 close the attacked source-side, attacked m2_psi4 target-side, and attacked direct m2 psi7 source-action blockers on the canonical-ontology-supported pre-observer lane only, and R47 reduces the remaining direct m2 psi7 source-eom role assignment witness to one still-missing coefficient-identification witness on fixed common support psi7(x), but that coefficient-identification witness remains absent, the target-slot assignment witness for m2_psi10 remains absent, the two remaining direct m2 pairwise witnesses remain absent, the direct g4, g6, and gY family zero witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_source_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "direct_m2_psi7_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "direct_m2_psi7_source_eom_common_monomial_support_packet_present": True,
            "direct_m2_psi7_source_eom_monomial_coefficient_identification_witness_present": False,
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
        "decomposition_of_R44_source_eom_gap": {
            "from_R44": "explicit_source_eom_role_assignment_witness_for_m2_psi7_to_mu_m2_plus3_segment_psi7_psi10",
            "reduced_by_R47_into": "explicit_source_eom_monomial_coefficient_identification_witness_for_m2_psi7_and_mu_m2_plus3_segment_psi7_psi10_on_common_psi7_of_x_support",
            "scope": "single_direct_m2_source_eom_role_common_local_support_only",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P60",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "decomposition_of_R44_source_eom_gap": artifact["decomposition_of_R44_source_eom_gap"],
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
