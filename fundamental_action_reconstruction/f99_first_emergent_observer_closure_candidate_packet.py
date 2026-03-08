#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f99_first_emergent_observer_closure_candidate_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f98 = load_json(
        "fundamental_action_reconstruction/generated/f98_first_emergent_observer_fixed_point_object_candidate_packet_summary.json"
    )

    closure_map = [[1.0]]
    fixed_object_vector = [float(x) for x in f98["source_observer_fixed_point_object"]["output_vector"]]
    closure_vector = mat_vec(closure_map, fixed_object_vector)

    summary = {
        "stage": "F99",
        "lane": "first_emergent_observer_closure_candidate_only",
        "goal": "export_one_emergent_observer_closure_candidate_without_claiming_actual_closure",
        "status": "F99_EXECUTED_FIRST_EMERGENT_OBSERVER_CLOSURE_CANDIDATE_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_fixed_point_object_input": "M_obs_fixed_object_preLM_v1",
        "observer_closure_candidate_map": {
            "object": "N_obs_closure_candidate_preLM_v1",
            "domain_basis": ["p_fix"],
            "codomain_basis": ["c_closure"],
            "matrix": closure_map,
            "closure_candidate_dimension": 1,
        },
        "observer_closure_candidate_properties": {
            "derived_only_from_fixed_point_object_state": True,
            "strict_core_only": True,
            "downstream_closure_candidate_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_closure_candidate": {
            "input_basis": ["p_fix"],
            "input_vector": fixed_object_vector,
            "output_basis": ["c_closure"],
            "output_vector": closure_vector,
            "c_closure_v1": closure_vector[0],
            "positive_closure_candidate_amplitude": closure_vector[0] > 0.0,
            "one_dimensional_closure_candidate": len(closure_vector) == 1,
        },
        "closure_candidate_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
