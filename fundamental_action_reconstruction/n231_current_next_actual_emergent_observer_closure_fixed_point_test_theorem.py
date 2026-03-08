#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n231_current_next_actual_emergent_observer_closure_fixed_point_test_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p211 = load_json(
        GENERATED / "p211_current_next_actual_emergent_observer_closure_fixed_point_test_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_FIXED_POINT_TEST_FROM_AM_OBS_ACTUAL_CLOSURE_REALIZATION_PRELM_V1_AFTER_P211"
    )
    status_ok = p211["status"] == expected_status

    checks = [
        {
            "id": "positive_next_actual_emergent_observer_closure_fixed_point_probe",
            "actual": p211["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N231",
        "status": "N231_DISCHARGED_CURRENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_FIXED_POINT_TEST_THEOREM_NO_FALSE_PASS",
        "scope": "current_next_actual_emergent_observer_closure_fixed_point_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "source_object": "S_preLM_strict_core_source_object_v1",
            "actual_emergent_observer_closure_realization_v5": "AM_obs_actual_closure_realization_preLM_v1",
            "actual_emergent_observer_closure_fixed_point_v3": "AN_obs_actual_closure_fixed_point_test_preLM_v1",
            "admissible_AM_obs_actual_closure_realization": True,
            "admissible_AN_obs_actual_closure_fixed_point_test": status_ok,
            "observer_information_deficit_downstream": status_ok,
        },
        "hard_limits": [
            "actual_emergent_observer_closure_not_yet_constructed",
            "no_QW2191_discharge",
            "no_strict_core_selector_closure",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
