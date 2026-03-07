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
    / "p43_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_source_eom_coherence_instance.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p43_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_source_eom_coherence_instance_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p42 = load_json(
        "fundamental_action_reconstruction/generated/p42_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_source_eom_coefficient_defect_polynomial_packet.json"
    )
    ax11 = load_json("fundamental_action_reconstruction/generated/ax11_canonical_ontology_supported_preobserver_source_eom_coherence_instance.json")

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
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
            "id": "p42_route_negative_before_ax11",
            "actual": p42["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_ONLY_ATTACKED_SOURCE_ACTION_BLOCKER_CLOSED_AND_SOURCE_EOM_DEFECT_POLYNOMIAL_EXPORTED_ROUTE_STILL_NOT_CLOSED_AFTER_R34",
            "meaning": "before AX11, the route still lacked the source-eom defect zero witness",
        },
        {
            "id": "ax11_local_source_eom_closure_present",
            "actual": ax11["result"]["canonical_ontology_supported_source_eom_coefficient_defect_zero_witness_present"],
            "expected": True,
            "meaning": "AX11 closes the attacked source-eom zero witness on the canonical-ontology-supported lane only",
        },
        {
            "id": "ax11_strict_core_promotion_forbidden",
            "actual": ax11["result"]["strict_core_promotion"],
            "expected": False,
            "meaning": "AX11 remains outside strict core",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "P43",
        "goal": "rerun_the_direct_formal_pair1_c1s1_family_route_after_closing_the_attacked_R34_source_eom_blocker_on_the_canonical_ontology_supported_preobserver_lane_only",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_ONLY_ATTACKED_SOURCE_ACTION_AND_SOURCE_EOM_BLOCKERS_CLOSED_ROUTE_STILL_NOT_CLOSED_AFTER_AX11",
        "reason": "AX10 closes the attacked direct m2 psi1 source-action blocker and AX11 closes the attacked direct m2 psi1 source-eom blocker on the canonical-ontology-supported pre-observer lane only, but the target-side assignment witness for m2_psi4 remains absent, the direct g4, g6, and gY family zero witnesses remain absent, the other direct m2 pairwise witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "lane": "canonical_ontology_supported_direct_formal_c1s1_family_route_after_AX10_and_AX11_only",
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_source_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "strict_core_source_side_blockers_closed": False,
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
        "remaining_missing_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P43",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
