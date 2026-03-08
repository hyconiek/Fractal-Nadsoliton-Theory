#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F128 = GENERATED / "f128_first_source_topology_selector_promotion_target_packet_summary.json"
IN_F146 = GENERATED / "f146_first_actual_source_topology_full_nontriviality_witness_packet_summary.json"
IN_F88 = GENERATED / "f88_first_exported_preobserver_source_object_orientation_export_packet_summary.json"
IN_F89 = GENERATED / "f89_first_exported_preobserver_selector_bridge_operator_packet_summary.json"
IN_F90 = GENERATED / "f90_first_exported_preobserver_selector_reduction_operator_packet_summary.json"
IN_F91 = GENERATED / "f91_first_exported_preobserver_selector_output_operator_packet_summary.json"
IN_N196 = GENERATED / "n196_current_exported_preobserver_source_object_admissible_orientation_export_theorem_summary.json"
IN_N197 = GENERATED / "n197_current_exported_preobserver_selector_bridge_operator_theorem_summary.json"
IN_N198 = GENERATED / "n198_current_exported_preobserver_selector_reduction_operator_theorem_summary.json"
IN_N199 = GENERATED / "n199_current_exported_preobserver_selector_output_operator_theorem_summary.json"
IN_N234 = GENERATED / "n234_current_global_selector_closure_and_qw2191_discharge_promotion_obstruction_theorem_summary.json"
OUT = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f128 = load_json(IN_F128)
    f146 = load_json(IN_F146)
    f88 = load_json(IN_F88)
    f89 = load_json(IN_F89)
    f90 = load_json(IN_F90)
    f91 = load_json(IN_F91)
    n196 = load_json(IN_N196)
    n197 = load_json(IN_N197)
    n198 = load_json(IN_N198)
    n199 = load_json(IN_N199)
    n234 = load_json(IN_N234)

    preobserver_scope_realized = (
        f90["selector_reduction_properties"]["preobserver_only"] is True
        and f91["selector_output_properties"]["preobserver_only"] is True
        and n234["theorem_result"]["observer_downstream_only"] is True
    )

    actual_selector_witness_exported = (
        f128["input_packet"] == "tau_src_candidate_v1"
        and f128["codomain_packet"]["name"] == "Sigma_sel_src_target_v1"
        and f146["input_packet"] == "tau_src_candidate_v1"
        and f146["full_source_topology_nontriviality_discharged"] is True
        and n196["theorem_result"]["admissible_E_orient"] is True
        and n197["theorem_result"]["admissible_B_sel"] is True
        and n198["theorem_result"]["admissible_R_sel"] is True
        and n199["theorem_result"]["admissible_O_sel"] is True
        and f89["source_alignment_witness"]["positive_signed_selector_response"] is True
        and f91["source_selector_output_response"]["positive_plus_output"] is True
        and f91["source_selector_output_response"]["vanishing_minus_output"] is True
        and preobserver_scope_realized is True
        and f128["observer_role"] == "downstream_only"
    )

    support_packet = {
        "input_packet": "tau_src_candidate_v1",
        "full_nontriviality_witness": f146["witness"],
        "selector_axis_realization": {
            "object": "E_orient_preLM_v1",
            "selector_axis_v1": f88["exported_orientation"]["selector_axis_v1"],
            "frame_basis": f88["exported_orientation"]["frame_basis"],
        },
        "selector_signed_split_realization": {
            "object": "B_sel_preLM_v1",
            "selector_projectors": f89["selector_bridge_operator"]["selector_projectors"],
            "positive_signed_selector_response": f89["source_alignment_witness"]["positive_signed_selector_response"],
        },
        "preobserver_scope_realization": {
            "selector_reduction_operator": "R_sel_preLM_v1",
            "selector_output_operator": "O_sel_preLM_v1",
            "preobserver_only": preobserver_scope_realized,
            "positive_plus_output": f91["source_selector_output_response"]["positive_plus_output"],
            "vanishing_minus_output": f91["source_selector_output_response"]["vanishing_minus_output"],
            "observer_downstream_only": n234["theorem_result"]["observer_downstream_only"],
        },
    }

    summary = {
        "packet_id": "F147",
        "status": "F147_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "refines_future_selector_target": f128["promotion_target"],
        "input_packet": "tau_src_candidate_v1",
        "witness": "Pi_sel_src_actual_witness_v1",
        "codomain_packet": f128["codomain_packet"]["name"],
        "support_packet_id": "W_src_selector_support_packet_v1",
        "support_packet": support_packet,
        "observer_role": "downstream_only",
        "chart_bound_preobserver_realization": True,
        "tau_src_identified_with_s_prelm": False,
        "actual_full_source_topology_nontriviality_witness_exported": f146["actual_full_source_topology_nontriviality_witness_exported"],
        "full_source_topology_nontriviality_discharged": f146["full_source_topology_nontriviality_discharged"],
        "actual_selector_witness_exported": actual_selector_witness_exported,
        "basis_independence_discharged": False,
        "qw2191_quotient_safe_discharged": False,
        "current_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "legacy_to_strict_bridge_claimed": False,
        "no_false_pass": True,
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
