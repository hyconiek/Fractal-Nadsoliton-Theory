#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f121_next_actual_emergent_observer_closure_commit_map_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f120 = load_json(
        GENERATED / "f120_next_actual_emergent_observer_closure_object_candidate_packet_summary.json"
    )

    x_actual_closure_obj_v2 = f120["source_observer_actual_closure_object_candidate"][
        "x_actual_closure_obj_v2"
    ]
    matrix = [[1.0], [0.0]]
    y_actual_commit_v2 = x_actual_closure_obj_v2
    y_actual_residual_v2 = 0.0

    summary = {
        "stage": "F121",
        "status": "F121_EXECUTED_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_COMMIT_MAP_PACKET_NO_FALSE_PASS",
        "input_object": "AK_obs_actual_closure_object_candidate_preLM_v1",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_commit_operator": {
            "operator": "AL_obs_actual_closure_commit_preLM_v1",
            "domain_basis": ["x_actual_closure_obj"],
            "codomain_basis": ["y_actual_commit", "y_actual_residual"],
            "matrix": matrix,
            "actual_closure_commit_dimension": 2,
        },
        "source_observer_actual_closure_commit": {
            "input_basis": ["x_actual_closure_obj"],
            "input_vector": [x_actual_closure_obj_v2],
            "output_basis": ["y_actual_commit", "y_actual_residual"],
            "output_vector": [y_actual_commit_v2, y_actual_residual_v2],
            "y_actual_commit_v2": y_actual_commit_v2,
            "y_actual_residual_v2": y_actual_residual_v2,
            "positive_actual_closure_commit_amplitude": y_actual_commit_v2 > 0.0,
            "actual_closure_residual_annihilated": y_actual_residual_v2 == 0.0,
        },
        "observer_actual_closure_commit_properties": {
            "derived_only_from_actual_closure_object_candidate_state": True,
            "strict_core_only": True,
            "downstream_actual_closure_commit_only": True,
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
