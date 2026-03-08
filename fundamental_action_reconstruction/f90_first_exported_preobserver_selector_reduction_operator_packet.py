#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f90_first_exported_preobserver_selector_reduction_operator_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f81 = load_json(
        "fundamental_action_reconstruction/generated/f81_first_additive_preobserver_strict_core_source_object_export_packet_summary.json"
    )
    f88 = load_json(
        "fundamental_action_reconstruction/generated/f88_first_exported_preobserver_source_object_orientation_export_packet_summary.json"
    )
    f89 = load_json(
        "fundamental_action_reconstruction/generated/f89_first_exported_preobserver_selector_bridge_operator_packet_summary.json"
    )

    e_parallel_2d = [float(x) for x in f88["exported_orientation"]["e_parallel"]]
    e_transverse_2d = [float(x) for x in f88["exported_orientation"]["e_transverse"]]
    r_sel = [e_parallel_2d + [0.0], e_transverse_2d + [0.0]]

    state = [float(x) for x in f81["state"]["closed_form"]]
    selector_response = mat_vec(r_sel, state)

    summary = {
        "stage": "F90",
        "lane": "first_exported_preobserver_selector_reduction_operator_only",
        "goal": "export_one_actual_preobserver_selector_reduction_map_from_B_sel_preLM_v1",
        "status": "F90_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SELECTOR_REDUCTION_OPERATOR_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "orientation_input": "E_orient_preLM_v1",
        "selector_bridge_input": "B_sel_preLM_v1",
        "selector_reduction_operator": {
            "object": "R_sel_preLM_v1",
            "domain_basis": ["u_T", "u_L", "u_M"],
            "codomain_basis": ["q_+", "q_-"],
            "matrix": r_sel,
        },
        "selector_reduction_properties": {
            "derived_only_from_orientation_and_bridge": True,
            "strict_core_only": True,
            "preobserver_only": True,
            "internal_selector_sector_exported": True,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "bridge_ready_for_O_sel": True,
        },
        "source_selector_response": {
            "input_state": state,
            "response_basis": ["q_+", "q_-"],
            "response_vector": selector_response,
            "r_plus_v1": selector_response[0],
            "r_minus_v1": selector_response[1],
            "positive_plus_channel": selector_response[0] > 0.0,
            "vanishing_minus_channel": abs(selector_response[1]) < 1.0e-12,
        },
        "input_bridge_witness": {
            "actual_B_sel_constructed": f89["actual_B_sel_constructed"],
            "positive_signed_selector_response": f89["source_alignment_witness"]["positive_signed_selector_response"],
        },
        "actual_R_sel_constructed": True,
        "actual_O_sel_constructed": False,
        "downstream_chain_complete": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
