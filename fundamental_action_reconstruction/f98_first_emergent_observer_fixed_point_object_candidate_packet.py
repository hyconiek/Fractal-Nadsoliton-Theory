#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f98_first_emergent_observer_fixed_point_object_candidate_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f97 = load_json(
        "fundamental_action_reconstruction/generated/f97_first_emergent_observer_fixed_point_reduction_packet_summary.json"
    )

    fixed_object_map = [[1.0]]
    fixed_point_vector = [float(x) for x in f97["source_observer_fixed_point"]["output_vector"]]
    fixed_object_vector = mat_vec(fixed_object_map, fixed_point_vector)

    summary = {
        "stage": "F98",
        "lane": "first_emergent_observer_fixed_point_object_candidate_only",
        "goal": "export_one_emergent_observer_fixed_point_object_candidate_without_claiming_actual_observer",
        "status": "F98_EXECUTED_FIRST_EMERGENT_OBSERVER_FIXED_POINT_OBJECT_CANDIDATE_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_fixed_point_input": "K_obs_fixed_point_preLM_v1",
        "observer_fixed_point_object_map": {
            "object": "M_obs_fixed_object_preLM_v1",
            "domain_basis": ["f_commit"],
            "codomain_basis": ["p_fix"],
            "matrix": fixed_object_map,
            "fixed_point_object_dimension": 1,
        },
        "observer_fixed_point_object_properties": {
            "derived_only_from_fixed_point_state": True,
            "strict_core_only": True,
            "downstream_fixed_point_object_only": True,
            "actual_emergent_observer_constructed": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_fixed_point_object": {
            "input_basis": ["f_commit"],
            "input_vector": fixed_point_vector,
            "output_basis": ["p_fix"],
            "output_vector": fixed_object_vector,
            "p_fix_v1": fixed_object_vector[0],
            "positive_fixed_point_object_amplitude": fixed_object_vector[0] > 0.0,
            "one_dimensional_fixed_point_object": len(fixed_object_vector) == 1,
        },
        "fixed_point_object_candidate_exported": True,
        "actual_emergent_observer_constructed": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
