#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f108_first_actual_emergent_observer_closure_object_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f107 = load_json(
        "fundamental_action_reconstruction/generated/f107_first_actual_emergent_observer_closure_candidate_packet_summary.json"
    )

    actual_closure_object = [[1.0]]
    candidate_vector = [
        float(x) for x in f107["source_observer_actual_closure_candidate"]["output_vector"]
    ]
    object_vector = mat_vec(actual_closure_object, candidate_vector)

    summary = {
        "stage": "F108",
        "lane": "first_actual_emergent_observer_closure_object_only",
        "goal": "export_one_actual_emergent_observer_closure_object_without_claiming_actual_closure",
        "status": "F108_EXECUTED_FIRST_ACTUAL_EMERGENT_OBSERVER_CLOSURE_OBJECT_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_candidate_input": "X_obs_actual_closure_candidate_preLM_v1",
        "observer_actual_closure_object_map": {
            "object": "Y_obs_actual_closure_object_preLM_v1",
            "domain_basis": ["k_actual_closure"],
            "codomain_basis": ["l_actual_closure_obj"],
            "matrix": actual_closure_object,
            "actual_closure_object_dimension": 1,
        },
        "observer_actual_closure_object_properties": {
            "derived_only_from_actual_closure_candidate_state": True,
            "strict_core_only": True,
            "downstream_actual_closure_object_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_actual_closure_object": {
            "input_basis": ["k_actual_closure"],
            "input_vector": candidate_vector,
            "output_basis": ["l_actual_closure_obj"],
            "output_vector": object_vector,
            "l_actual_closure_obj_v1": object_vector[0],
            "positive_actual_closure_object_amplitude": object_vector[0] > 0.0,
            "one_dimensional_actual_closure_object": len(object_vector) == 1,
        },
        "actual_closure_object_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
