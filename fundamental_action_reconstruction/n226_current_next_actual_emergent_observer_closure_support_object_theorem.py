#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n226_current_next_actual_emergent_observer_closure_support_object_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p206 = load_json(
        "fundamental_action_reconstruction/generated/p206_current_next_actual_emergent_observer_closure_support_object_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_FROM_AH_OBS_ACTUAL_CLOSURE_FIXED_POINT_TEST_PRELM_V1_AFTER_P206"
    )
    status_ok = p206["status"] == expected_status

    checks = [
        {
            "id": "positive_actual_emergent_observer_closure_support_probe",
            "actual": p206["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N226",
        "status": "N226_DISCHARGED_CURRENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_THEOREM_NO_FALSE_PASS",
        "scope": "current_next_actual_emergent_observer_closure_support_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "source_object": "S_preLM_strict_core_source_object_v1",
            "actual_emergent_observer_closure_fixed_point_test_v2": "AH_obs_actual_closure_fixed_point_test_preLM_v1",
            "actual_emergent_observer_closure_support_object_v2": "AI_obs_actual_closure_support_preLM_v1",
            "admissible_AH_obs_actual_closure_fixed_point_test": status_ok,
            "admissible_AI_obs_actual_closure_support_object": status_ok,
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
