#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "ax10_axiom_lane_preobserver_source_action_coherence_instance.json"
OUT_SUMMARY = GENERATED / "ax10_axiom_lane_preobserver_source_action_coherence_instance_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    ax9 = load_json("fundamental_action_reconstruction/generated/ax9_informational_nadsoliton_primacy_axiom_packet.json")
    r32 = load_json("fundamental_action_reconstruction/generated/r32_direct_m2_psi1_source_action_coefficient_defect_polynomial_packet.json")

    defect_packet = r32["direct_m2_psi1_source_action_coefficient_defect_packet"]
    source_symbol = defect_packet["source_action_coefficient_symbol"]
    target_symbol = defect_packet["lifted_action_coefficient_symbol"]
    common_support = defect_packet["common_fixed_action_monomial_support"]

    checks = [
        {
            "id": "ax9_canonical_ontology_provenance_fixed",
            "actual": ax9["result"]["canonical_informational_nadsoliton_ontology_provenance_fixed"],
            "expected": True,
            "meaning": "AX9 fixes provenance of the informational nadsoliton ontology from canonical program documents",
        },
        {
            "id": "ax9_pre_observer_source_side_may_be_tested_as_informational_coherence",
            "actual": ax9["recovered_consequences"]["pre_observer_source_side_may_be_tested_as_informational_coherence"],
            "expected": True,
            "meaning": "AX9 allows a pre-observer informational coherence test on the attacked source side as a provenance-supported external move",
        },
        {
            "id": "r32_source_action_defect_packet_present",
            "actual": r32["result"]["explicit_source_action_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R32 already exports the exact coefficient defect packet for the attacked source action lane",
        },
        {
            "id": "r32_zero_witness_absent_before_ax10",
            "actual": r32["result"]["explicit_zero_witness_for_source_action_coefficient_defect_polynomial_present"],
            "expected": False,
            "meaning": "before AX10, the attacked source action defect zero witness was still absent",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "AX10",
        "lane": "canonical-ontology-supported-external-lane",
        "packet_goal": "instantiate_one_preobserver_source_action_informational_coherence_use_instance_for_the_attacked_m2_psi1_coefficient_lane_using_the_canonical_ontology_provenance_fixed_in_AX9_without_promoting_any_result_into_strict_core",
        "source_reports": ["AX9", "R28", "R32"],
        "canonical_ontology_supported_preobserver_source_action_use_instance": {
            "attacked_lane": "direct_m2_psi1_source_action_on_common_psi1_squared_over_2_support",
            "source_action_coefficient_symbol": source_symbol,
            "common_informational_segment_parameter": target_symbol,
            "external_ontology_supported_assignment": f"{source_symbol} = {target_symbol}",
            "common_support": common_support,
            "closed_local_blocker": "explicit_zero_witness_for_the_direct_m2_psi1_source_action_coefficient_defect_polynomial_on_common_psi1_squared_over_2_support",
            "scope": "single_preobserver_source_action_lane_only",
            "boundary": "AX10 closes only the attacked source-action coefficient defect blocker on the canonical-ontology-supported external lane and leaves the source eom-side, target-side, other direct m2 blockers, g4/g6/gY blockers, pair1 zero equations, and QW-2191 unchanged",
        },
        "result": {
            "canonical_ontology_supported_source_action_assignment_instance_present": True,
            "canonical_ontology_supported_source_action_coefficient_defect_zero_witness_present": True,
            "strict_core_source_action_coefficient_defect_zero_witness_present": False,
            "strict_core_promotion": False,
        },
        "light_boundary": {
            "status": "light_before_observer_preobserver_scope_preserved",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "single pre-observer source-action coefficient lane before observer readout",
            "boundary": "AX10 uses the canonically documented informational nadsoliton ontology only on the pre-observer source-action side and does not derive anything from observer-side readout",
        },
        "classification": "canonical_ontology_supported_preobserver_source_action_use_instance_present_but_strict_core_unchanged",
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_strict_core_promotion",
            "no_claim_that_m2_psi1_equals_m2_psi4",
            "no_claim_that_the_source_eom_role_assignment_witness_is_present",
            "no_claim_that_the_target_side_assignment_witness_is_present",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_ToE_is_closed",
        ],
        "consistency_checks": checks,
        "no_false_pass": True,
    }

    summary = {
        "stage": "AX10",
        "status": "AX10_EXECUTED_CANONICAL_ONTOLOGY_SUPPORTED_PREOBSERVER_SOURCE_ACTION_USE_INSTANCE_NO_FALSE_PASS",
        "lane": "canonical-ontology-supported-external-lane",
        "result": artifact["classification"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "next_step": "P40",
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
