#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r9_existing_kernel_feedback_host_to_control_pushforward_packet_for_kobs.json"
OUT_SUMMARY = GENERATED / "r9_existing_kernel_feedback_host_to_control_pushforward_packet_for_kobs_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p12 = load_json(
        "fundamental_action_reconstruction/generated/p12_existing_kernel_feedback_factorization_rerun_after_host_carrier_packet.json"
    )
    r8 = load_json("fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs.json")
    r8_summary = load_json(
        "fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs_summary.json"
    )
    c14 = load_json("fundamental_action_reconstruction/generated/c14_control_mode_to_psi_transport_schema_summary.json")
    c15 = load_json("fundamental_action_reconstruction/generated/c15_control_only_pullback_submatrix_packet_summary.json")
    h8 = load_json("fundamental_action_reconstruction/generated/h8_minimal_component_carrier_construction_spec_summary.json")
    h33 = load_json("fundamental_action_reconstruction/generated/h33_pair1_selector_target_justification_audit.json")
    h34 = load_json("fundamental_action_reconstruction/generated/h34_basis_covariance_target_independence_audit.json")

    psi_basis = r8["declared_host_carrier"]["basis_order"]
    control_basis = c15["formal_objects"]["control_basis"]

    checks = [
        {
            "id": "p12_originally_missing_typed_projection",
            "actual": "typed_projection_or_pushforward_map_from_existing_kernel_feedback_into_the_explicit_H3_slot_chain_or_current_pair_block"
            in p12["remaining_missing_upstream_objects"],
            "expected": True,
            "meaning": "P12 indeed localized the typed projection/pushforward as one of the remaining factorization blockers",
        },
        {
            "id": "r8_host_carrier_present",
            "actual": r8_summary["result"],
            "expected": "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "meaning": "R8 already exports the host-side legacy carrier",
        },
        {
            "id": "c14_control_transport_schema_present",
            "actual": c14["result"]["control_transport_schema_present"],
            "expected": "yes",
            "meaning": "C14 already exports the control transport schema",
        },
        {
            "id": "c15_control_pullback_formula_present",
            "actual": c15["result"]["control_only_pullback_formula_present"],
            "expected": "yes",
            "meaning": "C15 already exports the control-only pullback formula built from that transport",
        },
        {
            "id": "c15_control_basis_matches_expected",
            "actual": control_basis,
            "expected": ["c1", "s1", "c2", "s2"],
            "meaning": "the control carrier basis is explicitly declared",
        },
        {
            "id": "h8_explicit_h3_chain_starts_on_pair1",
            "actual": h8["pair"]["plane"],
            "expected": "V_1 = span{c1, s1}",
            "meaning": "the explicit H3 chain still starts on pair1 rather than on the full control carrier",
        },
        {
            "id": "h33_pair1_still_local_chart_only",
            "actual": h33["result"],
            "expected": "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "meaning": "pair1 remains only a local chart and not a privileged selector target",
        },
        {
            "id": "h34_no_strict_covariance_argument",
            "actual": h34["status"],
            "expected": "PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT",
            "meaning": "no basis-covariant or target-independent selector reduction is present",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R9",
        "map_scope": "typed_host_to_control_pushforward_packet",
        "source_reports": ["R8", "C14", "C15", "H8", "H33", "H34", "P12"],
        "domain_carrier": {
            "label": "Psi_host_12",
            "basis_order": psi_basis,
            "dimension": len(psi_basis),
        },
        "codomain_carrier": {
            "label": "M_control",
            "basis_order": control_basis,
            "dimension": len(control_basis),
            "interpretation": "control-side mode carrier family preceding any selector-sector reduction",
        },
        "transport_schema": {
            "source_map_symbol": "T_control : B_control -> Psi_host_12",
            "source_map_shape": c15["formal_objects"]["T_control_shape"],
            "pushforward_map_symbol": "P_control = T_control^T : Psi_host_12 -> M_control",
            "pushforward_map_shape": [c15["formal_objects"]["T_control_shape"][1], c15["formal_objects"]["T_control_shape"][0]],
            "assembly_formula_on_control_side": c15["formal_objects"]["assembly_formula"],
        },
        "route_role": {
            "resolved_subobject": "typed_projection_or_pushforward_map_from_existing_kernel_feedback_into_the_explicit_H3_slot_chain_or_current_pair_block",
            "resolved_at_level": "host_to_control_carrier_pushforward_only",
            "unresolved_downstream_step": "selector_sector_reduction_from_M_control_to_pair1_or_an_equivalent_actual_target",
        },
        "formal_operator_pushforward_status": {
            "host_to_concrete_psi_block_identification_present": False,
            "reason": "C10_B1 remains open, so R9 exports the typed carrier map itself and not a coefficient-level evaluated pushforward of the host operator",
        },
        "reaches_explicit_current_pair_chain_directly": False,
        "selector_sector_reduction_present": False,
        "intertwiner_or_equality_witness_present": False,
        "control_scope_only": True,
        "consistency_checks": checks,
        "factorization_status": "typed_host_to_control_pushforward_present_but_no_legacy_selector_sector_reduction_to_pair1_or_equivalent_actual_target",
        "classification": "typed_host_to_control_pushforward_present_but_control_scope_only_and_not_yet_reduced_to_the_explicit_current_pair_chain",
        "frontier": "R9_B1",
        "no_false_pass": True,
    }

    summary = {
        "stage": "R9",
        "status": "PASS_PARTIAL_TYPED_HOST_TO_CONTROL_PUSHFORWARD_PACKET_READY_SELECTOR_REDUCTION_ABSENT",
        "result": artifact["classification"],
        "frontier": [
            "R9_B1",
            "H33_B1",
            "H34_B1",
            "H15_B1",
            "H16_B1",
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
