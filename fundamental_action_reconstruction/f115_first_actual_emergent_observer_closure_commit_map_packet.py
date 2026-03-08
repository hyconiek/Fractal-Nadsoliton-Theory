#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f115_first_actual_emergent_observer_closure_commit_map_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f114 = load_json(GENERATED / "f114_first_actual_emergent_observer_closure_object_candidate_packet_summary.json")

    actual_closure_commit = [[1.0], [0.0]]
    object_vector = [
        float(x) for x in f114["source_observer_actual_closure_object_candidate"]["output_vector"]
    ]
    commit_vector = mat_vec(actual_closure_commit, object_vector)

    summary = {
        "stage": "F115",
        "lane": "first_actual_emergent_observer_closure_commit_only",
        "goal": "export_one_actual_emergent_observer_closure_commit_map_without_claiming_actual_closure",
        "status": "F115_EXECUTED_FIRST_ACTUAL_EMERGENT_OBSERVER_CLOSURE_COMMIT_MAP_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_object_candidate_input": "AE_obs_actual_closure_object_candidate_preLM_v1",
        "observer_actual_closure_commit_map": {
            "object": "AF_obs_actual_closure_commit_preLM_v1",
            "domain_basis": ["r_actual_closure_obj"],
            "codomain_basis": ["s_actual_commit", "s_actual_residual"],
            "matrix": actual_closure_commit,
            "actual_closure_commit_dimension": 2,
        },
        "observer_actual_closure_commit_properties": {
            "derived_only_from_actual_closure_object_candidate_state": True,
            "strict_core_only": True,
            "downstream_actual_closure_commit_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_actual_closure_commit": {
            "input_basis": ["r_actual_closure_obj"],
            "input_vector": object_vector,
            "output_basis": ["s_actual_commit", "s_actual_residual"],
            "output_vector": commit_vector,
            "s_actual_commit_v1": commit_vector[0],
            "s_actual_residual_v1": commit_vector[1],
            "positive_actual_closure_commit_amplitude": commit_vector[0] > 0.0,
            "zero_actual_closure_residual": commit_vector[1] == 0.0,
        },
        "actual_closure_commit_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
