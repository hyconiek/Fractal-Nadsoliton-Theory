#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f102_first_emergent_observer_closure_support_object_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f101 = load_json(
        "fundamental_action_reconstruction/generated/f101_first_emergent_observer_closure_fixed_point_test_packet_summary.json"
    )

    closure_support_map = [[1.0]]
    closure_fix_vector = [
        float(x) for x in f101["source_observer_closure_fixed_point"]["output_vector"]
    ]
    closure_support_vector = mat_vec(closure_support_map, closure_fix_vector)

    summary = {
        "stage": "F102",
        "lane": "first_emergent_observer_closure_support_object_only",
        "goal": "export_one_emergent_observer_closure_support_object_without_claiming_actual_closure",
        "status": "F102_EXECUTED_FIRST_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_closure_fixed_point_input": "R_obs_closure_fixed_point_test_preLM_v1",
        "observer_closure_support_map": {
            "object": "S_obs_closure_support_preLM_v1",
            "domain_basis": ["e_closure_fix"],
            "codomain_basis": ["f_closure_support"],
            "matrix": closure_support_map,
            "closure_support_dimension": 1,
        },
        "observer_closure_support_properties": {
            "derived_only_from_closure_fixed_point_state": True,
            "strict_core_only": True,
            "downstream_closure_support_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_closure_support": {
            "input_basis": ["e_closure_fix"],
            "input_vector": closure_fix_vector,
            "output_basis": ["f_closure_support"],
            "output_vector": closure_support_vector,
            "f_closure_support_v1": closure_support_vector[0],
            "positive_closure_support_amplitude": closure_support_vector[0] > 0.0,
            "one_dimensional_closure_support": len(closure_support_vector) == 1,
        },
        "closure_support_object_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
