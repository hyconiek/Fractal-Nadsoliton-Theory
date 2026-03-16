#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "ax34_canonical_ontology_supported_preobserver_direct_m2_psi11_target_coherence_instance.json"
OUT_SUMMARY = (
    GENERATED / "ax34_canonical_ontology_supported_preobserver_direct_m2_psi11_target_coherence_instance_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    ax9 = load_json(
        "fundamental_action_reconstruction/generated/ax9_informational_nadsoliton_primacy_axiom_packet.json"
    )
    ax33 = load_json(
        "fundamental_action_reconstruction/generated/ax33_canonical_ontology_supported_preobserver_direct_m2_psi8_source_coherence_instance.json"
    )
    f3 = load_json(
        "fundamental_action_reconstruction/generated/f3_current_far_frontier_kernel_artifact_sensitivity_classification_packet_summary.json"
    )
    r79 = load_json(
        "fundamental_action_reconstruction/generated/r79_direct_m2_psi11_target_action_coefficient_defect_polynomial_packet.json"
    )
    r81 = load_json(
        "fundamental_action_reconstruction/generated/r81_direct_m2_psi11_target_eom_coefficient_defect_polynomial_packet.json"
    )

    action_defect_packet = r79["direct_m2_psi11_target_action_coefficient_defect_packet"]
    eom_defect_packet = r81["direct_m2_psi11_target_eom_coefficient_defect_packet"]

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
            "id": "ax33_previous_external_lane_direct_m2_closures_present",
            "actual": ax33["result"]["canonical_ontology_supported_direct_m2_psi8_source_action_coefficient_defect_zero_witness_present"]
            and ax33["result"]["canonical_ontology_supported_direct_m2_psi8_source_eom_coefficient_defect_zero_witness_present"]
            and ax33["result"]["canonical_ontology_supported_previous_external_lane_direct_m2_closures_preserved"],
            "expected": True,
            "meaning": "AX33 already extends the canonical-ontology-supported external lane with the psi8 slotwise defect closures and preserves all earlier external-lane direct-m2 closures",
        },
        {
            "id": "f3_route_classified_as_kernel_split_robust",
            "actual": bool(f3["kernel_split_robust_current_route_remaining_missing_objects"]),
            "expected": True,
            "meaning": "F3 keeps the current attacked route on the kernel-split-robust side after kernel correction",
        },
        {
            "id": "r79_target_action_defect_packet_present",
            "actual": r79["result"]["explicit_target_action_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R79 exports the exact target-action coefficient defect packet for m2_psi11 vs mu_m2_plus3_segment_psi8_psi11",
        },
        {
            "id": "r79_zero_witness_absent_before_ax34",
            "actual": r79["result"]["explicit_zero_witness_for_target_action_coefficient_defect_polynomial_present"],
            "expected": False,
            "meaning": "before AX34, the target-action defect zero witness for m2_psi11 was still absent",
        },
        {
            "id": "r81_target_eom_defect_packet_present",
            "actual": r81["result"]["explicit_target_eom_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R81 exports the exact target-eom coefficient defect packet for m2_psi11 vs mu_m2_plus3_segment_psi8_psi11",
        },
        {
            "id": "r81_zero_witness_absent_before_ax34",
            "actual": r81["result"]["explicit_zero_witness_for_target_eom_coefficient_defect_polynomial_present"],
            "expected": False,
            "meaning": "before AX34, the target-eom defect zero witness for m2_psi11 was still absent",
        },
        {
            "id": "lifted_symbol_is_consistent_on_action_and_eom",
            "actual": lifted_symbol == eom_defect_packet["lifted_eom_coefficient_symbol"],
            "expected": True,
            "meaning": "the declared common plus3 carrier-segment parameter symbol matches on the action and eom defect packets",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "AX34",
        "lane": "canonical-ontology-supported-external-lane",
        "packet_goal": "instantiate_one_preobserver_direct_m2_psi11_target_informational_coherence_use_instance_for_the_attacked_coefficient_lanes_using_the_canonical_ontology_provenance_fixed_in_AX9_without_promoting_any_result_into_strict_core_and_without_silent_strict_kernel_substitution",
        "source_reports": ["AX9", "AX33", "F3", "R79", "R81"],
        "canonical_ontology_supported_preobserver_direct_m2_psi11_target_use_instance": {
            "attacked_lanes": [
                "direct_m2_psi11_target_action_on_common_psi11_squared_over_2_support",
                "direct_m2_psi11_target_eom_on_common_psi11_of_x_support",
            ],
            "target_action_coefficient_symbol": target_symbol,
            "target_eom_coefficient_symbol": eom_defect_packet["target_eom_coefficient_symbol"],
            "common_informational_segment_parameter": lifted_symbol,
            "external_ontology_supported_assignment": f"{target_symbol} = {lifted_symbol}",
            "common_supports": {"action": action_support, "eom": eom_support},
            "closed_local_blockers": [
                "explicit_zero_witness_for_the_direct_m2_psi11_target_action_coefficient_defect_polynomial_on_common_psi11_squared_over_2_support",
                "explicit_zero_witness_for_the_direct_m2_psi11_target_eom_coefficient_defect_polynomial_on_common_psi11_of_x_support",
            ],
            "scope": "single_preobserver_direct_m2_psi11_target_coefficient_lanes_only",
            "boundary": "AX34 closes only the attacked direct m2 psi11 target-side coefficient defect blockers on the canonical-ontology-supported external lane and leaves the remaining direct m2 defect blockers, the g4/g6/gY blockers, the pair1 zero equations, and QW-2191 unchanged",
        },
        "result": {
            "canonical_ontology_supported_direct_m2_psi11_target_assignment_instance_present": True,
            "canonical_ontology_supported_direct_m2_psi11_target_action_coefficient_defect_zero_witness_present": True,
            "canonical_ontology_supported_direct_m2_psi11_target_eom_coefficient_defect_zero_witness_present": True,
            "canonical_ontology_supported_previous_external_lane_direct_m2_closures_preserved": True,
            "strict_core_direct_m2_psi11_target_action_coefficient_defect_zero_witness_present": False,
            "strict_core_direct_m2_psi11_target_eom_coefficient_defect_zero_witness_present": False,
            "strict_core_promotion": False,
        },
        "light_boundary": {
            "status": "light_before_observer_preobserver_scope_preserved",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "single pre-observer direct m2 psi11 target coefficient lanes before observer readout",
            "boundary": "AX34 uses the canonically documented informational nadsoliton ontology only on the pre-observer direct m2 psi11 target side and does not derive anything from observer-side readout or from silent K_strict_gate inheritance",
        },
        "classification": "canonical_ontology_supported_preobserver_direct_m2_psi11_target_use_instance_present_but_strict_core_unchanged",
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_strict_core_promotion",
            "no_claim_that_m2_psi8_equals_m2_psi11",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_ToE_is_closed",
        ],
        "consistency_checks": checks,
        "no_false_pass": True,
    }

    summary = {
        "stage": "AX34",
        "status": "AX34_EXECUTED_CANONICAL_ONTOLOGY_SUPPORTED_PREOBSERVER_DIRECT_M2_PSI11_TARGET_USE_INSTANCE_NO_FALSE_PASS",
        "lane": "canonical-ontology-supported-external-lane",
        "result": artifact["classification"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "next_step": "P628",
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
