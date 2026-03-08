#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f122_next_actual_emergent_observer_closure_realization_map_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f121 = load_json(
        GENERATED / "f121_next_actual_emergent_observer_closure_commit_map_packet_summary.json"
    )

    y_actual_commit_v2 = f121["source_observer_actual_closure_commit"]["y_actual_commit_v2"]
    y_actual_residual_v2 = f121["source_observer_actual_closure_commit"]["y_actual_residual_v2"]
    matrix = [[1.0, 0.0]]
    z_actual_closure_real_v2 = y_actual_commit_v2

    summary = {
        "stage": "F122",
        "status": "F122_EXECUTED_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_REALIZATION_MAP_PACKET_NO_FALSE_PASS",
        "input_object": "AL_obs_actual_closure_commit_preLM_v1",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_realization_operator": {
            "operator": "AM_obs_actual_closure_realization_preLM_v1",
            "domain_basis": ["y_actual_commit", "y_actual_residual"],
            "codomain_basis": ["z_actual_closure_real"],
            "matrix": matrix,
            "actual_closure_realization_dimension": 1,
        },
        "source_observer_actual_closure_realization": {
            "input_basis": ["y_actual_commit", "y_actual_residual"],
            "input_vector": [y_actual_commit_v2, y_actual_residual_v2],
            "output_basis": ["z_actual_closure_real"],
            "output_vector": [z_actual_closure_real_v2],
            "z_actual_closure_real_v2": z_actual_closure_real_v2,
            "positive_actual_closure_realization_amplitude": z_actual_closure_real_v2 > 0.0,
            "actual_closure_residual_annihilated": y_actual_residual_v2 == 0.0,
        },
        "observer_actual_closure_realization_properties": {
            "derived_only_from_actual_closure_commit_state": True,
            "strict_core_only": True,
            "downstream_actual_closure_realization_only": True,
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
