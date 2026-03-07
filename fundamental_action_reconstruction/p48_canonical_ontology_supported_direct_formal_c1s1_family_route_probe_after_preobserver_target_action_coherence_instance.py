#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p48_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_target_action_coherence_instance.json"
OUT_SUMMARY = GENERATED / "p48_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_target_action_coherence_instance_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p46 = load_json(
        "fundamental_action_reconstruction/generated/p46_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_action_coefficient_defect_polynomial_packet.json"
    )
    ax12 = load_json(
        "fundamental_action_reconstruction/generated/ax12_canonical_ontology_supported_preobserver_target_action_coherence_instance.json"
    )

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
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
            "id": "p46_target_action_defect_witness_was_still_missing",
            "actual": "explicit_zero_witness_for_the_direct_m2_psi4_target_action_coefficient_defect_polynomial_on_common_psi4_squared_over_2_support"
            in p46["remaining_missing_objects"],
            "expected": True,
            "meaning": "P46 still had the attacked target-action defect zero witness as a missing object",
        },
        {
            "id": "ax12_target_action_defect_witness_closed_on_external_lane",
            "actual": ax12["result"]["canonical_ontology_supported_target_action_coefficient_defect_zero_witness_present"],
            "expected": True,
            "meaning": "AX12 closes the attacked target-action defect zero witness on the canonical-ontology-supported lane",
        },
        {
            "id": "ax12_strict_core_target_action_witness_still_absent",
            "actual": ax12["result"]["strict_core_target_action_coefficient_defect_zero_witness_present"],
            "expected": False,
            "meaning": "AX12 does not promote the attacked target-action closure into strict core",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "P48",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-AX11-and-AX12",
        "goal": "rerun_the_canonical_ontology_supported_direct_family_route_after_local_closure_of_the_attacked_m2_psi4_target_action_defect_witness",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_AND_ATTACKED_TARGET_ACTION_BLOCKERS_CLOSED_ROUTE_STILL_NOT_CLOSED_AFTER_AX12",
        "reason": "AX10 and AX11 locally close the attacked source-side blockers on the canonical-ontology-supported pre-observer lane, AX12 locally closes the attacked m2_psi4 target-action coefficient-defect blocker on that same lane, but the target eom-role assignment witness remains absent, the direct g4, g6, and gY family zero witnesses remain absent, the other direct m2 pairwise witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_source_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
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
        "remaining_missing_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P48",
        "status": artifact["status"],
        "lane": artifact["lane"],
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
