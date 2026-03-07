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
    / "p51_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_target_eom_coherence_instance.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p51_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_target_eom_coherence_instance_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p50 = load_json(
        "fundamental_action_reconstruction/generated/p50_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet.json"
    )
    ax13 = load_json(
        "fundamental_action_reconstruction/generated/ax13_canonical_ontology_supported_preobserver_target_eom_coherence_instance.json"
    )

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10",
        "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
        "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    checks = [
        {
            "id": "p50_target_eom_defect_polynomial_was_still_missing_before_ax13",
            "actual": p50["route_state"]["target_eom_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX13, P50 still recorded the target-eom defect zero witness as missing",
        },
        {
            "id": "ax13_local_target_eom_closure_present",
            "actual": ax13["result"]["canonical_ontology_supported_target_eom_coefficient_defect_zero_witness_present"],
            "expected": True,
            "meaning": "AX13 adds the local closure of the attacked target-eom blocker on the canonical-ontology-supported lane",
        },
        {
            "id": "ax13_target_action_closure_preserved",
            "actual": ax13["result"]["canonical_ontology_supported_target_action_coefficient_defect_zero_witness_preserved"],
            "expected": True,
            "meaning": "AX13 preserves the already local target-action closure from AX12",
        },
        {
            "id": "ax13_source_side_closures_preserved",
            "actual": ax13["result"]["canonical_ontology_supported_source_eom_coefficient_defect_zero_witness_preserved"],
            "expected": True,
            "meaning": "AX13 preserves the already local source-side closures from AX10/AX11",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "P51",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-AX11-AX12-and-AX13",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_local_target_eom_closure_of_the_attacked_m2_psi4_blocker",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_AND_ATTACKED_TARGET_SIDE_BLOCKERS_CLOSED_ROUTE_STILL_NOT_CLOSED_AFTER_AX13",
        "reason": "AX10 and AX11 close the attacked source-side blockers on the canonical-ontology-supported pre-observer lane only, AX12 closes the attacked m2_psi4 target-action blocker on that same lane, AX13 closes the attacked m2_psi4 target-eom blocker on that same lane, but the direct g4, g6, and gY family zero witnesses remain absent, the remaining direct m2 pairwise witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_source_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
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
        "closed_local_blockers": [
            "R32_B1",
            "R34_B1",
            "R37_B1",
            "R39_B1",
        ],
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P51",
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
