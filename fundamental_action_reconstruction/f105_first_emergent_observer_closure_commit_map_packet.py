#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f105_first_emergent_observer_closure_commit_map_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f104 = load_json(
        "fundamental_action_reconstruction/generated/f104_first_emergent_observer_closure_object_candidate_packet_summary.json"
    )

    closure_commit = [
        [1.0],
        [0.0],
    ]
    object_vector = [float(x) for x in f104["source_observer_closure_object_candidate"]["output_vector"]]
    commit_vector = mat_vec(closure_commit, object_vector)

    summary = {
        "stage": "F105",
        "lane": "first_emergent_observer_closure_commit_only",
        "goal": "export_one_emergent_observer_closure_commit_map_without_claiming_actual_closure",
        "status": "F105_EXECUTED_FIRST_EMERGENT_OBSERVER_CLOSURE_COMMIT_MAP_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_closure_object_candidate_input": "U_obs_closure_object_candidate_preLM_v1",
        "observer_closure_commit_map": {
            "object": "V_obs_closure_commit_preLM_v1",
            "domain_basis": ["h_closure_obj"],
            "codomain_basis": ["i_commit", "i_residual"],
            "matrix": closure_commit,
            "closure_commit_dimension": 2,
        },
        "observer_closure_commit_properties": {
            "derived_only_from_closure_object_candidate_state": True,
            "strict_core_only": True,
            "downstream_closure_commit_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_closure_commit": {
            "input_basis": ["h_closure_obj"],
            "input_vector": object_vector,
            "output_basis": ["i_commit", "i_residual"],
            "output_vector": commit_vector,
            "i_commit_v1": commit_vector[0],
            "i_residual_v1": commit_vector[1],
            "positive_commit_amplitude": commit_vector[0] > 0.0,
            "zero_residual_channel": abs(commit_vector[1]) < 1e-12,
        },
        "closure_commit_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
