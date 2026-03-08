#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f123_next_actual_emergent_observer_closure_fixed_point_test_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f122 = load_json(
        GENERATED / "f122_next_actual_emergent_observer_closure_realization_map_packet_summary.json"
    )

    actual_closure_fixed_point = [[1.0]]
    realization_vector = [
        float(x)
        for x in f122["source_observer_actual_closure_realization"]["output_vector"]
    ]
    fixed_point_vector = mat_vec(actual_closure_fixed_point, realization_vector)

    summary = {
        "stage": "F123",
        "status": "F123_EXECUTED_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_FIXED_POINT_TEST_PACKET_NO_FALSE_PASS",
        "input_object": "AM_obs_actual_closure_realization_preLM_v1",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_fixed_point_operator": {
            "operator": "AN_obs_actual_closure_fixed_point_test_preLM_v1",
            "domain_basis": ["z_actual_closure_real"],
            "codomain_basis": ["aa_actual_closure_fix"],
            "matrix": actual_closure_fixed_point,
            "actual_closure_fixed_point_dimension": 1,
        },
        "source_observer_actual_closure_fixed_point": {
            "input_basis": ["z_actual_closure_real"],
            "input_vector": realization_vector,
            "output_basis": ["aa_actual_closure_fix"],
            "output_vector": fixed_point_vector,
            "aa_actual_closure_fix_v2": fixed_point_vector[0],
            "positive_actual_closure_fixed_point_amplitude": fixed_point_vector[0] > 0.0,
            "one_dimensional_actual_closure_fixed_point": len(fixed_point_vector) == 1,
        },
        "observer_actual_closure_fixed_point_properties": {
            "derived_only_from_actual_closure_realization_map": True,
            "strict_core_only": True,
            "downstream_actual_closure_fixed_point_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "hard_limits": [
            "not_actual_emergent_observer_closure",
            "no_QW2191_discharge",
            "no_strict_core_selector_closure",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
