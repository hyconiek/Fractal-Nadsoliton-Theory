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
    / "p59_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi7_source_action_coherence_instance.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p59_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi7_source_action_coherence_instance_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p58 = load_json(
        "fundamental_action_reconstruction/generated/p58_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi7_source_action_coefficient_defect_polynomial_packet.json"
    )
    ax14 = load_json(
        "fundamental_action_reconstruction/generated/ax14_canonical_ontology_supported_preobserver_direct_m2_psi7_source_action_coherence_instance.json"
    )

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_source_eom_role_assignment_witness_for_m2_psi7_to_mu_m2_plus3_segment_psi7_psi10",
        "explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10",
        "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
        "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    checks = [
        {
            "id": "p58_source_action_defect_polynomial_was_still_missing_before_ax14",
            "actual": p58["route_state"]["direct_m2_psi7_source_action_coefficient_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX14, P58 still recorded the direct m2 psi7 source-action coefficient defect zero witness as missing",
        },
        {
            "id": "ax14_local_direct_m2_psi7_source_action_closure_present",
            "actual": ax14["result"]["canonical_ontology_supported_direct_m2_psi7_source_action_coefficient_defect_zero_witness_present"],
            "expected": True,
            "meaning": "AX14 adds the local closure of the attacked direct m2 psi7 source-action blocker on the canonical-ontology-supported lane",
        },
        {
            "id": "ax14_previous_target_side_closures_preserved",
            "actual": ax14["result"]["canonical_ontology_supported_target_eom_coefficient_defect_zero_witness_preserved"],
            "expected": True,
            "meaning": "AX14 preserves the already local attacked target-side closures from AX12 and AX13",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "P59",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-AX11-AX12-AX13-R40-R41-R42-R43-R44-R45-R46-and-AX14",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_local_source_action_closure_of_the_attacked_direct_m2_psi7_blocker",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_DIRECT_M2_PSI7_SOURCE_ACTION_BLOCKER_CLOSED_ROUTE_STILL_NOT_CLOSED_AFTER_AX14",
        "reason": "AX10, AX11, AX12, and AX13 close the attacked source-side and attacked m2_psi4 target-side blockers on the canonical-ontology-supported pre-observer lane only, R40-R46 reduce the direct pairwise target m2_psi7 = m2_psi10 to one attacked source-action defect zero-witness gap plus the remaining source-eom, target-slot, and other-route gaps, and AX14 closes only that one attacked direct m2 psi7 source-action defect blocker on the same external lane, but the source-eom role witness remains absent, the target-slot assignment witness for m2_psi10 remains absent, the two remaining direct m2 pairwise witnesses remain absent, the direct g4, g6, and gY family zero witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_source_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "direct_m2_psi7_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "direct_m2_psi7_source_eom_role_assignment_witness_present": False,
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
        "closed_local_blockers": [
            "R32_B1",
            "R34_B1",
            "R37_B1",
            "R39_B1",
            "R46_B1",
        ],
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P59",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "closed_local_blockers": artifact["closed_local_blockers"],
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
