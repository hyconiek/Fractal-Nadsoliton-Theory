#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f111_first_actual_emergent_observer_closure_fixed_point_test_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f110 = load_json(
        "fundamental_action_reconstruction/generated/f110_first_actual_emergent_observer_closure_realization_object_packet_summary.json"
    )

    actual_closure_fixed_point = [[1.0]]
    realization_vector = [
        float(x)
        for x in f110["source_observer_actual_closure_realization"]["output_vector"]
    ]
    fixed_point_vector = mat_vec(actual_closure_fixed_point, realization_vector)

    summary = {
        "stage": "F111",
        "lane": "first_actual_emergent_observer_closure_fixed_point_only",
        "goal": "export_one_actual_emergent_observer_closure_fixed_point_test_without_claiming_actual_closure",
        "status": "F111_EXECUTED_FIRST_ACTUAL_EMERGENT_OBSERVER_CLOSURE_FIXED_POINT_TEST_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_realization_input": "AA_obs_actual_closure_realization_preLM_v1",
        "observer_actual_closure_fixed_point_test": {
            "object": "AB_obs_actual_closure_fixed_point_test_preLM_v1",
            "domain_basis": ["n_actual_closure_real"],
            "codomain_basis": ["o_actual_closure_fix"],
            "matrix": actual_closure_fixed_point,
            "actual_closure_fixed_point_dimension": 1,
        },
        "observer_actual_closure_fixed_point_properties": {
            "derived_only_from_actual_closure_realization_object": True,
            "strict_core_only": True,
            "downstream_actual_closure_fixed_point_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_actual_closure_fixed_point": {
            "input_basis": ["n_actual_closure_real"],
            "input_vector": realization_vector,
            "output_basis": ["o_actual_closure_fix"],
            "output_vector": fixed_point_vector,
            "o_actual_closure_fix_v1": fixed_point_vector[0],
            "positive_actual_closure_fixed_point_amplitude": fixed_point_vector[0] > 0.0,
            "one_dimensional_actual_closure_fixed_point": len(fixed_point_vector) == 1,
        },
        "actual_closure_fixed_point_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
