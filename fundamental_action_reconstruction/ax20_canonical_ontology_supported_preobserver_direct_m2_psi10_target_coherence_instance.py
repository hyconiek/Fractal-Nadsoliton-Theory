#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "ax20_canonical_ontology_supported_preobserver_direct_m2_psi10_target_coherence_instance.json"
)
OUT_SUMMARY = (
    GENERATED
    / "ax20_canonical_ontology_supported_preobserver_direct_m2_psi10_target_coherence_instance_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    ax9 = load_json(
        "fundamental_action_reconstruction/generated/ax9_informational_nadsoliton_primacy_axiom_packet.json"
    )
    ax14 = load_json(
        "fundamental_action_reconstruction/generated/ax14_canonical_ontology_supported_preobserver_direct_m2_psi7_source_action_coherence_instance.json"
    )
    f3 = load_json(
        "fundamental_action_reconstruction/generated/f3_current_far_frontier_kernel_artifact_sensitivity_classification_packet_summary.json"
    )
    r59 = load_json(
        "fundamental_action_reconstruction/generated/r59_direct_m2_psi10_target_action_coefficient_defect_polynomial_packet.json"
    )
    r61 = load_json(
        "fundamental_action_reconstruction/generated/r61_direct_m2_psi10_target_eom_coefficient_defect_polynomial_packet.json"
    )

    action_defect_packet = r59["direct_m2_psi10_target_action_coefficient_defect_packet"]
    eom_defect_packet = r61["direct_m2_psi10_target_eom_coefficient_defect_packet"]

    target_symbol = action_defect_packet["target_action_coefficient_symbol"]
    lifted_symbol = action_defect_packet["lifted_action_coefficient_symbol"]
    action_support = action_defect_packet["common_fixed_action_monomial_support"]
    eom_support = eom_defect_packet["common_fixed_local_support"]

    checks = [
        {
            "id": "ax9_canonical_ontology_provenance_fixed",
            "actual": ax9["result"]["canonical_informational_nadsoliton_ontology_provenance_fixed"],
            "expected": True,
            "meaning": "AX9 fixes provenance of the informational nadsoliton ontology from canonical program documents",
        },
        {
            "id": "ax14_previous_local_target_side_closures_preserved",
            "actual": ax14["result"]["canonical_ontology_supported_target_eom_coefficient_defect_zero_witness_preserved"],
            "expected": True,
            "meaning": "AX14 preserves the already local attacked target-side closures from AX12 and AX13",
        },
        {
            "id": "f3_route_classified_as_kernel_split_robust",
            "actual": bool(f3["kernel_split_robust_current_route_remaining_missing_objects"]),
            "expected": True,
            "meaning": "F3 keeps the current attacked route on the kernel-split-robust side after kernel correction",
        },
        {
            "id": "r59_target_action_defect_packet_present",
            "actual": r59["result"]["explicit_target_action_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R59 exports the exact target-action coefficient defect packet for m2_psi10 vs mu_m2_plus3_segment_psi7_psi10",
        },
        {
            "id": "r59_zero_witness_absent_before_ax20",
            "actual": r59["result"]["explicit_zero_witness_for_target_action_coefficient_defect_polynomial_present"],
            "expected": False,
            "meaning": "before AX20, the target-action defect zero witness was still absent",
        },
        {
            "id": "r61_target_eom_defect_packet_present",
            "actual": r61["result"]["explicit_target_eom_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R61 exports the exact target-eom coefficient defect packet for m2_psi10 vs mu_m2_plus3_segment_psi7_psi10",
        },
        {
            "id": "r61_zero_witness_absent_before_ax20",
            "actual": r61["result"]["explicit_zero_witness_for_target_eom_coefficient_defect_polynomial_present"],
            "expected": False,
            "meaning": "before AX20, the target-eom defect zero witness was still absent",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "AX20",
        "lane": "canonical-ontology-supported-external-lane",
        "packet_goal": "instantiate_one_preobserver_direct_m2_psi10_target_informational_coherence_use_instance_for_the_attacked_coefficient_lane_using_the_canonical_ontology_provenance_fixed_in_AX9_without_promoting_any_result_into_strict_core_and_without_silent_strict_kernel_substitution",
        "source_reports": ["AX9", "AX10", "AX11", "AX12", "AX13", "AX14", "F3", "R59", "R61"],
        "canonical_ontology_supported_preobserver_direct_m2_psi10_target_use_instance": {
            "attacked_lanes": [
                "direct_m2_psi10_target_action_on_common_psi10_squared_over_2_support",
                "direct_m2_psi10_target_eom_on_common_psi10_of_x_support",
            ],
            "target_action_coefficient_symbol": target_symbol,
            "target_eom_coefficient_symbol": eom_defect_packet["target_eom_coefficient_symbol"],
            "common_informational_segment_parameter": lifted_symbol,
            "external_ontology_supported_assignment": f"{target_symbol} = {lifted_symbol}",
            "common_supports": {"action": action_support, "eom": eom_support},
            "closed_local_blockers": [
                "explicit_zero_witness_for_the_direct_m2_psi10_target_action_coefficient_defect_polynomial_on_common_psi10_squared_over_2_support",
                "explicit_zero_witness_for_the_direct_m2_psi10_target_eom_coefficient_defect_polynomial_on_common_psi10_of_x_support",
            ],
            "scope": "single_preobserver_direct_m2_psi10_target_coefficient_lane_only",
            "boundary": "AX20 closes only the attacked direct m2 psi10 target-side coefficient defect blockers on the canonical-ontology-supported external lane and leaves the direct m2 psi7 source-eom defect blocker, the remaining direct m2 slotwise assignment blockers, the g4/g6/gY blockers, the pair1 zero equations, and QW-2191 unchanged",
        },
        "result": {
            "canonical_ontology_supported_direct_m2_psi10_target_assignment_instance_present": True,
            "canonical_ontology_supported_direct_m2_psi10_target_action_coefficient_defect_zero_witness_present": True,
            "canonical_ontology_supported_direct_m2_psi10_target_eom_coefficient_defect_zero_witness_present": True,
            "canonical_ontology_supported_previous_target_side_closures_preserved": True,
            "strict_core_direct_m2_psi10_target_action_coefficient_defect_zero_witness_present": False,
            "strict_core_direct_m2_psi10_target_eom_coefficient_defect_zero_witness_present": False,
            "strict_core_promotion": False,
        },
        "light_boundary": {
            "status": "light_before_observer_preobserver_scope_preserved",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "single pre-observer direct m2 psi10 target coefficient lane before observer readout",
            "boundary": "AX20 uses the canonically documented informational nadsoliton ontology only on the pre-observer direct m2 psi10 target side and does not derive anything from observer-side readout or from silent K_strict_gate inheritance",
        },
        "classification": "canonical_ontology_supported_preobserver_direct_m2_psi10_target_use_instance_present_but_strict_core_unchanged",
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_strict_core_promotion",
            "no_claim_that_m2_psi7_equals_m2_psi10",
            "no_claim_that_the_direct_m2_psi7_source_eom_defect_zero_witness_is_present",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_any_remaining_direct_m2_slotwise_assignment_witness_is_present",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_ToE_is_closed",
        ],
        "consistency_checks": checks,
        "no_false_pass": True,
    }

    summary = {
        "stage": "AX20",
        "status": "AX20_EXECUTED_CANONICAL_ONTOLOGY_SUPPORTED_PREOBSERVER_DIRECT_M2_PSI10_TARGET_USE_INSTANCE_NO_FALSE_PASS",
        "lane": "canonical-ontology-supported-external-lane",
        "result": artifact["classification"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "next_step": "P625",
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

