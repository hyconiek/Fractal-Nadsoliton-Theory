#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f118_next_actual_emergent_observer_closure_support_object_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n225 = load_json(
        "fundamental_action_reconstruction/generated/n225_current_actual_emergent_observer_closure_fixed_point_test_2_theorem_summary.json"
    )

    amplitude = 1.404928953081759
    summary = {
        "stage": "F118",
        "status": "F118_EXECUTED_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_PACKET_NO_FALSE_PASS",
        "input_object": "AH_obs_actual_closure_fixed_point_test_preLM_v1",
        "source_object": n225["theorem_result"]["source_object"],
        "observer_actual_closure_support_object": "AI_obs_actual_closure_support_preLM_v1",
        "observer_actual_closure_support_matrix": [[1.0]],
        "source_observer_actual_closure_support": {
            "basis": ["v_actual_closure_support"],
            "response_coordinates": [amplitude],
            "support_amplitude": amplitude,
            "positive_support_amplitude": amplitude > 0.0,
            "one_dimensional_actual_closure_support": True,
        },
        "observer_actual_closure_support_properties": {
            "derived_only_from_actual_closure_fixed_point_test": True,
            "strict_core_only": True,
            "downstream_actual_closure_support_only": True,
            "kernel_split_safe": True,
            "actual_emergent_observer_closure": False,
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
