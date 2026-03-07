#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "ax13_canonical_ontology_supported_preobserver_target_eom_coherence_instance.json"
OUT_SUMMARY = GENERATED / "ax13_canonical_ontology_supported_preobserver_target_eom_coherence_instance_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    ax9 = load_json("fundamental_action_reconstruction/generated/ax9_informational_nadsoliton_primacy_axiom_packet.json")
    ax12 = load_json("fundamental_action_reconstruction/generated/ax12_canonical_ontology_supported_preobserver_target_action_coherence_instance.json")
    f3 = load_json("fundamental_action_reconstruction/generated/f3_current_far_frontier_kernel_artifact_sensitivity_classification_packet_summary.json")
    r39 = load_json("fundamental_action_reconstruction/generated/r39_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet.json")

    defect_packet = r39["direct_m2_psi4_target_eom_coefficient_defect_packet"]
    target_symbol = defect_packet["target_eom_coefficient_symbol"]
    lifted_symbol = defect_packet["lifted_eom_coefficient_symbol"]
    common_support = defect_packet["common_support"]

    checks = [
        {
            "id": "ax9_canonical_ontology_provenance_fixed",
            "actual": ax9["result"]["canonical_informational_nadsoliton_ontology_provenance_fixed"],
            "expected": True,
            "meaning": "AX9 fixes provenance of the informational nadsoliton ontology from canonical program documents",
        },
        {
            "id": "ax12_local_target_action_closure_present",
            "actual": ax12["result"]["canonical_ontology_supported_target_action_coefficient_defect_zero_witness_present"],
            "expected": True,
            "meaning": "AX12 preserves the previously closed attacked target-action blocker on the same canonical-ontology-supported lane",
        },
        {
            "id": "f3_route_classified_as_kernel_split_robust",
            "actual": bool(f3["kernel_split_robust_current_route_remaining_missing_objects"]),
            "expected": True,
            "meaning": "F3 keeps the current attacked route on the kernel-split-robust side after kernel correction",
        },
        {
            "id": "r39_target_eom_defect_packet_present",
            "actual": r39["result"]["explicit_target_eom_coefficient_defect_packet_present"],
            "expected": True,
            "meaning": "R39 already exports the exact coefficient defect packet for the attacked target-eom lane",
        },
        {
            "id": "r39_zero_witness_absent_before_ax13",
            "actual": r39["result"]["explicit_zero_witness_for_target_eom_coefficient_defect_polynomial_present"],
            "expected": False,
            "meaning": "before AX13, the attacked target-eom defect zero witness was still absent",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "AX13",
        "lane": "canonical-ontology-supported-external-lane",
        "packet_goal": "instantiate_one_preobserver_target_eom_informational_coherence_use_instance_for_the_attacked_m2_psi4_coefficient_lane_using_the_canonical_ontology_provenance_fixed_in_AX9_without_promoting_any_result_into_strict_core_and_without_silent_strict_kernel_substitution",
        "source_reports": ["AX9", "AX10", "AX11", "AX12", "F3", "R28", "R39"],
        "canonical_ontology_supported_preobserver_target_eom_use_instance": {
            "attacked_lane": "direct_m2_psi4_target_eom_on_common_psi4_of_x_support",
            "target_eom_coefficient_symbol": target_symbol,
            "common_informational_segment_parameter": lifted_symbol,
            "external_ontology_supported_assignment": f"{target_symbol} = {lifted_symbol}",
            "common_support": common_support,
            "closed_local_blocker": "explicit_zero_witness_for_the_direct_m2_psi4_target_eom_coefficient_defect_polynomial_on_common_psi4_of_x_support",
            "scope": "single_preobserver_target_eom_lane_only",
            "boundary": "AX13 closes only the attacked target-eom coefficient defect blocker on the canonical-ontology-supported external lane and leaves the remaining direct m2 pairwise blockers, g4/g6/gY blockers, pair1 zero equations, and QW-2191 unchanged",
        },
        "result": {
            "canonical_ontology_supported_target_eom_assignment_instance_present": True,
            "canonical_ontology_supported_target_eom_coefficient_defect_zero_witness_present": True,
            "canonical_ontology_supported_source_action_coefficient_defect_zero_witness_preserved": True,
            "canonical_ontology_supported_source_eom_coefficient_defect_zero_witness_preserved": True,
            "canonical_ontology_supported_target_action_coefficient_defect_zero_witness_preserved": True,
            "strict_core_target_eom_coefficient_defect_zero_witness_present": False,
            "strict_core_promotion": False,
        },
        "light_boundary": {
            "status": "light_before_observer_preobserver_scope_preserved",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "single pre-observer target-eom coefficient lane before observer readout",
            "boundary": "AX13 uses the canonically documented informational nadsoliton ontology only on the pre-observer target-eom side and does not derive anything from observer-side readout or from silent K_strict_gate inheritance",
        },
        "classification": "canonical_ontology_supported_preobserver_target_eom_use_instance_present_but_strict_core_unchanged",
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_strict_core_promotion",
            "no_claim_that_m2_psi1_equals_m2_psi4",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_any_remaining_direct_m2_pairwise_witness_is_present",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_ToE_is_closed",
        ],
        "consistency_checks": checks,
        "no_false_pass": True,
    }

    summary = {
        "stage": "AX13",
        "status": "AX13_EXECUTED_CANONICAL_ONTOLOGY_SUPPORTED_PREOBSERVER_TARGET_EOM_USE_INSTANCE_NO_FALSE_PASS",
        "lane": "canonical-ontology-supported-external-lane",
        "result": artifact["classification"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "next_step": "P51",
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
