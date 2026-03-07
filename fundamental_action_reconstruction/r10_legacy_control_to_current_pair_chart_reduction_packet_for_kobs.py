#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r10_legacy_control_to_current_pair_chart_reduction_packet_for_kobs.json"
OUT_SUMMARY = GENERATED / "r10_legacy_control_to_current_pair_chart_reduction_packet_for_kobs_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p13 = load_json(
        "fundamental_action_reconstruction/generated/p13_existing_kernel_feedback_factorization_rerun_after_host_to_control_pushforward_packet.json"
    )
    r9 = load_json("fundamental_action_reconstruction/generated/r9_existing_kernel_feedback_host_to_control_pushforward_packet_for_kobs.json")
    r9_summary = load_json(
        "fundamental_action_reconstruction/generated/r9_existing_kernel_feedback_host_to_control_pushforward_packet_for_kobs_summary.json"
    )
    c15 = load_json("fundamental_action_reconstruction/generated/c15_control_only_pullback_submatrix_packet_summary.json")
    c29 = load_json("fundamental_action_reconstruction/generated/c29_local_projector_formula_packet_summary.json")
    h8 = load_json("fundamental_action_reconstruction/generated/h8_minimal_component_carrier_construction_spec_summary.json")
    h33 = load_json("fundamental_action_reconstruction/generated/h33_pair1_selector_target_justification_audit.json")
    h34 = load_json("fundamental_action_reconstruction/generated/h34_basis_covariance_target_independence_audit.json")

    control_basis = r9["codomain_carrier"]["basis_order"]
    pair_basis = h8["pair"]["basis"]
    pi_pair1 = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
    ]

    checks = [
        {
            "id": "p13_originally_missing_legacy_selector_sector_reduction",
            "actual": "selector_sector_reduction_of_the_legacy_control_side_onto_pair1_or_an_equivalent_actual_target"
            in p13["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P13 indeed localized the legacy selector-sector reduction as one remaining factorization blocker",
        },
        {
            "id": "r9_control_carrier_present",
            "actual": r9_summary["result"],
            "expected": "typed_host_to_control_pushforward_present_but_control_scope_only_and_not_yet_reduced_to_the_explicit_current_pair_chain",
            "meaning": "R9 already exports the legacy control-side carrier family",
        },
        {
            "id": "c15_control_basis_matches_expected",
            "actual": c15["formal_objects"]["control_basis"],
            "expected": ["c1", "s1", "c2", "s2"],
            "meaning": "the control basis is explicitly declared as (c1,s1,c2,s2)",
        },
        {
            "id": "h8_pair_basis_matches_expected",
            "actual": pair_basis,
            "expected": ["c1", "s1"],
            "meaning": "the explicit H3 chain starts on pair1 basis (c1,s1)",
        },
        {
            "id": "control_basis_prefix_matches_pair_basis",
            "actual": control_basis[:2],
            "expected": pair_basis,
            "meaning": "the chosen pair1 basis matches the first control-basis block and therefore admits a typed chart extraction",
        },
        {
            "id": "c29_local_projector_formula_present",
            "actual": c29["result"]["explicit_local_projector_formula_present"],
            "expected": "yes_partial",
            "meaning": "C29 already exports an explicit local reduced-projector formula on each pair plane",
        },
        {
            "id": "h33_pair1_local_chart_only",
            "actual": h33["result"],
            "expected": "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "meaning": "pair1 remains only a local chart and not a privileged selector target",
        },
        {
            "id": "h34_no_covariance_argument",
            "actual": h34["status"],
            "expected": "PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT",
            "meaning": "no basis-covariant or target-independent selector reduction is present",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R10",
        "map_scope": "legacy_control_to_current_pair_chart_reduction_packet",
        "source_reports": ["R9", "C15", "C29", "H8", "H33", "H34", "P13"],
        "domain_carrier": {
            "label": "M_control",
            "basis_order": control_basis,
            "dimension": len(control_basis),
        },
        "codomain_carrier": {
            "label": "V_1",
            "basis_order": pair_basis,
            "plane": h8["pair"]["plane"],
            "dimension": len(pair_basis),
        },
        "chart_reduction_map": {
            "symbol": "Pi_pair1 : M_control -> V_1",
            "matrix": pi_pair1,
            "construction_rule": "coordinate_block_extraction_onto_the_first_declared_control_pair_matching_the_explicit_current_pair_chain_domain",
        },
        "supporting_local_reduced_projector_family": {
            "source": "C29",
            "P_red(theta)": c29["serialized_formulas"]["P_red(theta)"],
            "status": "local_pair_level_only",
        },
        "route_role": {
            "resolved_subobject": "selector_sector_reduction_of_the_legacy_control_side_onto_pair1_or_an_equivalent_actual_target",
            "resolved_at_level": "chosen_explicit_current_pair_chart_only",
            "unresolved_downstream_step": "intertwiner_or_equality_witness_identifying_the_chart_reduced_legacy_side_with_the_computed_current_pair_H3_block",
        },
        "strict_selector_target_justification_present": False,
        "basis_covariant_target_independent_reduction_present": False,
        "intertwiner_or_equality_witness_present": False,
        "chart_scope_only": True,
        "consistency_checks": checks,
        "factorization_status": "legacy_control_to_chosen_current_pair_chart_reduction_present_but_no_intertwiner_or_equality_witness_to_the_computed_current_pair_block",
        "classification": "explicit_current_pair_chart_reduction_present_but_chart_scoped_only_and_not_a_strict_selector_target_justification",
        "frontier": "R10_B1",
        "no_false_pass": True,
    }

    summary = {
        "stage": "R10",
        "status": "PASS_PARTIAL_EXPLICIT_CURRENT_PAIR_CHART_REDUCTION_PACKET_READY_EQUALITY_WITNESS_ABSENT",
        "result": artifact["classification"],
        "frontier": [
            "R10_B1",
            "H16_B1",
            "H33_B1",
            "H34_B1",
            "C10_B1",
        ],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
