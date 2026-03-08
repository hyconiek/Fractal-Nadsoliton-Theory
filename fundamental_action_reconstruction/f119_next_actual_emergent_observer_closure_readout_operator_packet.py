#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f119_next_actual_emergent_observer_closure_readout_operator_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f118 = load_json(
        "fundamental_action_reconstruction/generated/f118_next_actual_emergent_observer_closure_support_object_packet_summary.json"
    )

    v_actual_closure_support_v1 = f118["source_observer_actual_closure_support"]["support_amplitude"]
    matrix = [[1.0], [0.0]]
    w_actual_commit_v1 = v_actual_closure_support_v1
    w_actual_gap_v1 = 0.0

    summary = {
        "stage": "F119",
        "status": "F119_EXECUTED_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_PACKET_NO_FALSE_PASS",
        "input_object": "AI_obs_actual_closure_support_preLM_v1",
        "source_object": f118["source_object"],
        "observer_actual_closure_readout_operator": {
            "operator": "AJ_obs_actual_closure_readout_preLM_v1",
            "domain_basis": ["v_actual_closure_support"],
            "codomain_basis": ["w_actual_commit", "w_actual_gap"],
            "matrix": matrix,
            "actual_closure_readout_dimension": 2,
        },
        "source_observer_actual_closure_readout": {
            "input_basis": ["v_actual_closure_support"],
            "input_vector": [v_actual_closure_support_v1],
            "output_basis": ["w_actual_commit", "w_actual_gap"],
            "output_vector": [w_actual_commit_v1, w_actual_gap_v1],
            "w_actual_commit_v1": w_actual_commit_v1,
            "w_actual_gap_v1": w_actual_gap_v1,
            "positive_commit_amplitude": w_actual_commit_v1 > 0.0,
            "zero_gap_channel": w_actual_gap_v1 == 0.0,
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
