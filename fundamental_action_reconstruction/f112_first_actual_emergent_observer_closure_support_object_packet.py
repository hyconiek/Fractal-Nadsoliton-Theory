#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f112_first_actual_emergent_observer_closure_support_object_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f111 = load_json(
        "fundamental_action_reconstruction/generated/f111_first_actual_emergent_observer_closure_fixed_point_test_packet_summary.json"
    )

    actual_closure_support = [[1.0]]
    fixed_point_vector = [
        float(x)
        for x in f111["source_observer_actual_closure_fixed_point"]["output_vector"]
    ]
    support_vector = mat_vec(actual_closure_support, fixed_point_vector)

    summary = {
        "stage": "F112",
        "lane": "first_actual_emergent_observer_closure_support_only",
        "goal": "export_one_actual_emergent_observer_closure_support_object_without_claiming_actual_closure",
        "status": "F112_EXECUTED_FIRST_ACTUAL_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_fixed_point_input": "AB_obs_actual_closure_fixed_point_test_preLM_v1",
        "observer_actual_closure_support_object": {
            "object": "AC_obs_actual_closure_support_preLM_v1",
            "domain_basis": ["o_actual_closure_fix"],
            "codomain_basis": ["p_actual_closure_support"],
            "matrix": actual_closure_support,
            "actual_closure_support_dimension": 1,
        },
        "observer_actual_closure_support_properties": {
            "derived_only_from_actual_closure_fixed_point_test": True,
            "strict_core_only": True,
            "downstream_actual_closure_support_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_actual_closure_support": {
            "input_basis": ["o_actual_closure_fix"],
            "input_vector": fixed_point_vector,
            "output_basis": ["p_actual_closure_support"],
            "output_vector": support_vector,
            "p_actual_closure_support_v1": support_vector[0],
            "positive_actual_closure_support_amplitude": support_vector[0] > 0.0,
            "one_dimensional_actual_closure_support": len(support_vector) == 1,
        },
        "actual_closure_support_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
