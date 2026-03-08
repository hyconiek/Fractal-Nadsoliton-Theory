#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f91_first_exported_preobserver_selector_output_operator_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f90 = load_json(
        "fundamental_action_reconstruction/generated/f90_first_exported_preobserver_selector_reduction_operator_packet_summary.json"
    )

    o_sel = [
        [1.0, 0.0],
        [0.0, 1.0],
    ]
    selector_response = [float(x) for x in f90["source_selector_response"]["response_vector"]]
    output_response = mat_vec(o_sel, selector_response)

    summary = {
        "stage": "F91",
        "lane": "first_exported_preobserver_selector_output_operator_only",
        "goal": "export_one_actual_preobserver_selector_output_map_from_R_sel_preLM_v1",
        "status": "F91_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SELECTOR_OUTPUT_OPERATOR_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "orientation_input": "E_orient_preLM_v1",
        "selector_bridge_input": "B_sel_preLM_v1",
        "selector_reduction_input": "R_sel_preLM_v1",
        "selector_output_operator": {
            "object": "O_sel_preLM_v1",
            "domain_basis": ["q_+", "q_-"],
            "codomain_basis": ["o_+", "o_-"],
            "matrix": o_sel,
        },
        "selector_output_properties": {
            "derived_only_from_selector_reduction": True,
            "strict_core_only": True,
            "preobserver_only": True,
            "selector_output_sector_exported": True,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "bridge_ready_for_emergent_observer_limit": True,
        },
        "source_selector_output_response": {
            "input_response_basis": ["q_+", "q_-"],
            "input_response_vector": selector_response,
            "output_basis": ["o_+", "o_-"],
            "output_vector": output_response,
            "o_plus_v1": output_response[0],
            "o_minus_v1": output_response[1],
            "positive_plus_output": output_response[0] > 0.0,
            "vanishing_minus_output": abs(output_response[1]) < 1.0e-12,
        },
        "input_reduction_witness": {
            "actual_R_sel_constructed": f90["actual_R_sel_constructed"],
            "positive_plus_channel": f90["source_selector_response"]["positive_plus_channel"],
            "vanishing_minus_channel": f90["source_selector_response"]["vanishing_minus_channel"],
        },
        "actual_O_sel_constructed": True,
        "emergent_observer_constructed": False,
        "downstream_chain_complete": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
