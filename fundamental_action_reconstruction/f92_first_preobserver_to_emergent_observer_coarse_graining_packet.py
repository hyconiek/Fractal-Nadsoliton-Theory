#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f92_first_preobserver_to_emergent_observer_coarse_graining_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f91 = load_json(
        "fundamental_action_reconstruction/generated/f91_first_exported_preobserver_selector_output_operator_packet_summary.json"
    )

    coarse = [
        [0.5, -0.5],
        [0.5, 0.5],
    ]
    output_vector = [float(x) for x in f91["source_selector_output_response"]["output_vector"]]
    coarse_vector = mat_vec(coarse, output_vector)

    summary = {
        "stage": "F92",
        "lane": "first_preobserver_to_emergent_observer_coarse_graining_only",
        "goal": "export_one_downstream_coarse_graining_map_from_Q_out_v1_without_claiming_actual_observer",
        "status": "F92_EXECUTED_FIRST_PREOBSERVER_TO_EMERGENT_OBSERVER_COARSE_GRAINING_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "selector_output_input": "O_sel_preLM_v1",
        "coarse_graining_operator": {
            "object": "C_obs_limit_preLM_v1",
            "domain_basis": ["o_+", "o_-"],
            "codomain_basis": ["y_bias", "y_total"],
            "matrix": coarse,
        },
        "coarse_graining_properties": {
            "derived_only_from_selector_output": True,
            "strict_core_only": True,
            "preobserver_to_observer_limit_only": True,
            "actual_emergent_observer_constructed": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_coarse_response": {
            "input_basis": ["o_+", "o_-"],
            "input_vector": output_vector,
            "output_basis": ["y_bias", "y_total"],
            "output_vector": coarse_vector,
            "y_bias_v1": coarse_vector[0],
            "y_total_v1": coarse_vector[1],
            "positive_bias": coarse_vector[0] > 0.0,
            "positive_total": coarse_vector[1] > 0.0,
        },
        "actual_coarse_graining_constructed": True,
        "emergent_observer_constructed": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
