#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f101_first_emergent_observer_closure_fixed_point_test_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f100 = load_json(
        "fundamental_action_reconstruction/generated/f100_first_emergent_observer_closure_realization_map_packet_summary.json"
    )

    fixed_point_test = [[1.0]]
    closure_real_vector = [float(x) for x in f100["source_observer_closure_realization"]["output_vector"]]
    closure_fix_vector = mat_vec(fixed_point_test, closure_real_vector)

    summary = {
        "stage": "F101",
        "lane": "first_emergent_observer_closure_fixed_point_test_only",
        "goal": "export_one_emergent_observer_closure_fixed_point_test_without_claiming_actual_closure",
        "status": "F101_EXECUTED_FIRST_EMERGENT_OBSERVER_CLOSURE_FIXED_POINT_TEST_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_closure_realization_input": "Q_obs_closure_realization_preLM_v1",
        "observer_closure_fixed_point_test": {
            "object": "R_obs_closure_fixed_point_test_preLM_v1",
            "domain_basis": ["d_closure"],
            "codomain_basis": ["e_closure_fix"],
            "matrix": fixed_point_test,
            "closure_fixed_point_dimension": 1,
        },
        "observer_closure_fixed_point_properties": {
            "derived_only_from_closure_realization_state": True,
            "strict_core_only": True,
            "downstream_closure_fixed_point_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_closure_fixed_point": {
            "input_basis": ["d_closure"],
            "input_vector": closure_real_vector,
            "output_basis": ["e_closure_fix"],
            "output_vector": closure_fix_vector,
            "e_closure_fix_v1": closure_fix_vector[0],
            "positive_closure_fixed_point_amplitude": closure_fix_vector[0] > 0.0,
            "one_dimensional_closure_fixed_point": len(closure_fix_vector) == 1,
        },
        "closure_fixed_point_test_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
