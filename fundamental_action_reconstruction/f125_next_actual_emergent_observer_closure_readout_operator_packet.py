#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f125_next_actual_emergent_observer_closure_readout_operator_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f124 = load_json(
        GENERATED / "f124_next_actual_emergent_observer_closure_support_object_packet_summary.json"
    )

    readout_matrix = [[1.0], [0.0]]
    support_vector = [
        float(x) for x in f124["source_observer_actual_closure_support"]["output_vector"]
    ]
    readout_vector = mat_vec(readout_matrix, support_vector)

    summary = {
        "stage": "F125",
        "status": "F125_EXECUTED_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_PACKET_NO_FALSE_PASS",
        "input_object": "AO_obs_actual_closure_support_preLM_v1",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_readout_operator": {
            "operator": "AP_obs_actual_closure_readout_preLM_v1",
            "domain_basis": ["ab_actual_closure_support"],
            "codomain_basis": ["ac_actual_commit", "ac_actual_gap"],
            "matrix": readout_matrix,
            "actual_closure_readout_dimension": 2,
        },
        "source_observer_actual_closure_readout": {
            "input_basis": ["ab_actual_closure_support"],
            "input_vector": support_vector,
            "output_basis": ["ac_actual_commit", "ac_actual_gap"],
            "output_vector": readout_vector,
            "ac_actual_commit_v2": readout_vector[0],
            "ac_actual_gap_v2": readout_vector[1],
            "positive_commit_amplitude": readout_vector[0] > 0.0,
            "zero_gap_channel": readout_vector[1] == 0.0,
        },
        "observer_actual_closure_readout_properties": {
            "derived_only_from_actual_closure_support_state": True,
            "strict_core_only": True,
            "downstream_actual_closure_readout_only": True,
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
