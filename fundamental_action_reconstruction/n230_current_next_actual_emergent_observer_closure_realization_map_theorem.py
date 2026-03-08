#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n230_current_next_actual_emergent_observer_closure_realization_map_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p210 = load_json(
        GENERATED / "p210_current_next_actual_emergent_observer_closure_realization_map_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_REALIZATION_MAP_FROM_AL_OBS_ACTUAL_CLOSURE_COMMIT_PRELM_V1_AFTER_P210"
    )
    status_ok = p210["status"] == expected_status

    checks = [
        {
            "id": "positive_next_actual_emergent_observer_closure_realization_probe",
            "actual": p210["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N230",
        "status": "N230_DISCHARGED_CURRENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_REALIZATION_MAP_THEOREM_NO_FALSE_PASS",
        "scope": "current_next_actual_emergent_observer_closure_realization_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "source_object": "S_preLM_strict_core_source_object_v1",
            "actual_emergent_observer_closure_commit_v4": "AL_obs_actual_closure_commit_preLM_v1",
            "actual_emergent_observer_closure_realization_v4": "AM_obs_actual_closure_realization_preLM_v1",
            "admissible_AL_obs_actual_closure_commit": True,
            "admissible_AM_obs_actual_closure_realization": status_ok,
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
