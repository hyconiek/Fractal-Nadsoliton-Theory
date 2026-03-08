#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n232_current_next_actual_emergent_observer_closure_support_object_theorem_summary.json"

EXPECTED_STATUS = (
    "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_"
    "SUPPORT_OBJECT_FROM_AN_OBS_ACTUAL_CLOSURE_FIXED_POINT_TEST_PRELM_V1_AFTER_P212"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p212 = load_json(
        GENERATED / "p212_current_next_actual_emergent_observer_closure_support_object_probe_summary.json"
    )

    status_ok = p212["status"] == EXPECTED_STATUS
    summary = {
        "stage": "N232",
        "status": (
            "N232_DISCHARGED_CURRENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_THEOREM_NO_FALSE_PASS"
            if status_ok
            else "N232_BLOCKED_CURRENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_THEOREM"
        ),
        "theorem_result": {
            "source_object": "S_preLM_strict_core_source_object_v1",
            "actual_emergent_observer_closure_fixed_point_v3": "AN_obs_actual_closure_fixed_point_test_preLM_v1",
            "actual_emergent_observer_closure_support_v3": "AO_obs_actual_closure_support_preLM_v1",
            "admissible_AN_obs_actual_closure_fixed_point_test": status_ok,
            "admissible_AO_obs_actual_closure_support_object": status_ok,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "observer_information_deficit_downstream": status_ok,
            "kernel_split_safe": status_ok,
        },
        "depends_on": [
            "N231_CURRENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_FIXED_POINT_TEST_THEOREM",
            "P212_CURRENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_PROBE",
        ],
        "no_false_pass": status_ok,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
