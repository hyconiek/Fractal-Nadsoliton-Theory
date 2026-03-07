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
    / "p50_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p50_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p49 = load_json(
        "fundamental_action_reconstruction/generated/p49_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_eom_common_local_support_packet.json"
    )
    r39 = load_json("fundamental_action_reconstruction/generated/r39_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet.json")

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_zero_witness_for_the_direct_m2_psi4_target_eom_coefficient_defect_polynomial_on_common_psi4_of_x_support",
        "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10",
        "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
        "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    checks = [
        {
            "id": "p49_local_target_action_closure_present",
            "actual": p49["route_state"]["attacked_target_action_blocker_closed_on_canonical_ontology_supported_lane_only"],
            "expected": True,
            "meaning": "P49 still records the local closure of the attacked target-action blocker on the canonical-ontology-supported lane only",
        },
        {
            "id": "p49_target_eom_gap_reduced_before_r39",
            "actual": not p49["route_state"]["target_eom_monomial_coefficient_identification_witness_present"],
            "expected": True,
            "meaning": "before R39, the target eom-side gap was already reduced to one common-support coefficient-identification witness",
        },
        {
            "id": "r39_target_eom_defect_polynomial_present",
            "actual": r39["result"]["explicit_target_eom_coefficient_defect_packet_present"],
            "expected": True,
            "meaning": "R39 adds the exact coefficient-defect packet for the attacked target eom-side witness",
        },
        {
            "id": "r39_target_eom_zero_witness_still_absent",
            "actual": r39["result"]["explicit_zero_witness_for_target_eom_coefficient_defect_polynomial_present"],
            "expected": False,
            "meaning": "R39 does not claim that the target eom-side defect zero witness already holds",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "P50",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-AX11-AX12-R38-and-R39",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_reducing_the_attacked_target_eom_side_to_one_exact_defect_polynomial_zero_witness_gap",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_AND_ATTACKED_TARGET_ACTION_BLOCKERS_CLOSED_AND_TARGET_EOM_DEFECT_POLYNOMIAL_EXPORTED_ROUTE_STILL_NOT_CLOSED_AFTER_R39",
        "reason": "AX10 and AX11 close the attacked source-side blockers on the canonical-ontology-supported pre-observer lane only, AX12 closes the attacked m2_psi4 target-action blocker on that same lane, R38 reduces the attacked target eom-side witness to one coefficient-identification witness on common psi4(x) support, and R39 reduces that witness further to one exact defect-polynomial zero-witness gap, but the zero witness remains absent, the direct g4, g6, and gY family zero witnesses remain absent, the other direct m2 pairwise witnesses remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open",
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_source_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "target_eom_defect_polynomial_exported": True,
            "target_eom_defect_zero_witness_present": False,
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
        "decomposition_of_P49_target_eom_gap": {
            "from_P49": "explicit_target_eom_monomial_coefficient_identification_witness_for_m2_psi4_and_mu_m2_plus3_segment_psi1_psi4_on_common_psi4_of_x_support",
            "reduced_by_R39_into": "explicit_zero_witness_for_the_direct_m2_psi4_target_eom_coefficient_defect_polynomial_on_common_psi4_of_x_support",
            "scope": "single_direct_m2_target_eom_coefficient_gap_only",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P50",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "decomposition_of_P49_target_eom_gap": artifact["decomposition_of_P49_target_eom_gap"],
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
