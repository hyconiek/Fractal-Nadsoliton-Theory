#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f120_next_actual_emergent_observer_closure_object_candidate_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f119 = load_json(
        GENERATED / "f119_next_actual_emergent_observer_closure_readout_operator_packet_summary.json"
    )

    w_actual_commit_v1 = f119["source_observer_actual_closure_readout"]["w_actual_commit_v1"]
    w_actual_gap_v1 = f119["source_observer_actual_closure_readout"]["w_actual_gap_v1"]
    matrix = [[1.0, 0.0]]
    x_actual_closure_obj_v2 = w_actual_commit_v1

    summary = {
        "stage": "F120",
        "status": "F120_EXECUTED_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_OBJECT_CANDIDATE_PACKET_NO_FALSE_PASS",
        "input_object": "AJ_obs_actual_closure_readout_preLM_v1",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_object_candidate_operator": {
            "operator": "AK_obs_actual_closure_object_candidate_preLM_v1",
            "domain_basis": ["w_actual_commit", "w_actual_gap"],
            "codomain_basis": ["x_actual_closure_obj"],
            "matrix": matrix,
            "actual_closure_object_candidate_dimension": 1,
        },
        "source_observer_actual_closure_object_candidate": {
            "input_basis": ["w_actual_commit", "w_actual_gap"],
            "input_vector": [w_actual_commit_v1, w_actual_gap_v1],
            "output_basis": ["x_actual_closure_obj"],
            "output_vector": [x_actual_closure_obj_v2],
            "x_actual_closure_obj_v2": x_actual_closure_obj_v2,
            "positive_actual_closure_object_candidate_amplitude": x_actual_closure_obj_v2 > 0.0,
            "actual_gap_annihilated": w_actual_gap_v1 == 0.0,
        },
        "observer_actual_closure_object_candidate_properties": {
            "derived_only_from_actual_closure_readout_state": True,
            "strict_core_only": True,
            "downstream_actual_closure_object_candidate_only": True,
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
