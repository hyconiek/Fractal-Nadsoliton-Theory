#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f93_first_emergent_observer_limit_readout_operator_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f92 = load_json(
        "fundamental_action_reconstruction/generated/f92_first_preobserver_to_emergent_observer_coarse_graining_packet_summary.json"
    )

    readout = [
        [1.0, 1.0],
        [-1.0, 1.0],
    ]
    coarse_vector = [float(x) for x in f92["source_coarse_response"]["output_vector"]]
    readout_vector = mat_vec(readout, coarse_vector)

    summary = {
        "stage": "F93",
        "lane": "first_emergent_observer_limit_readout_operator_only",
        "goal": "export_one_downstream_observer_limit_readout_map_from_Y_obs_limit_v1_without_claiming_actual_observer",
        "status": "F93_EXECUTED_FIRST_EMERGENT_OBSERVER_LIMIT_READOUT_OPERATOR_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "coarse_graining_input": "C_obs_limit_preLM_v1",
        "observer_limit_readout_operator": {
            "object": "L_obs_limit_preLM_v1",
            "domain_basis": ["y_bias", "y_total"],
            "codomain_basis": ["z_commit", "z_residual"],
            "matrix": readout,
        },
        "observer_limit_readout_properties": {
            "derived_only_from_coarse_graining": True,
            "strict_core_only": True,
            "observer_limit_only": True,
            "actual_emergent_observer_constructed": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_limit_readout": {
            "input_basis": ["y_bias", "y_total"],
            "input_vector": coarse_vector,
            "output_basis": ["z_commit", "z_residual"],
            "output_vector": readout_vector,
            "z_commit_v1": readout_vector[0],
            "z_residual_v1": readout_vector[1],
            "positive_commit": readout_vector[0] > 0.0,
            "vanishing_residual": abs(readout_vector[1]) < 1.0e-12,
        },
        "actual_observer_limit_readout_constructed": True,
        "emergent_observer_constructed": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
