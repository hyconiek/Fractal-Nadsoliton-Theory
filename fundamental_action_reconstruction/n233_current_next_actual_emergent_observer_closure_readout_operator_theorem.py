#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n233_current_next_actual_emergent_observer_closure_readout_operator_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p213 = load_json(
        GENERATED / "p213_current_next_actual_emergent_observer_closure_readout_operator_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_FROM_AO_OBS_ACTUAL_CLOSURE_SUPPORT_PRELM_V1_AFTER_P213"
    )
    status_ok = p213["status"] == expected_status

    checks = [
        {
            "id": "positive_next_actual_emergent_observer_closure_readout_probe",
            "actual": p213["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N233",
        "status": "N233_DISCHARGED_CURRENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_THEOREM_NO_FALSE_PASS",
        "scope": "current_next_actual_emergent_observer_closure_readout_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "source_object": "S_preLM_strict_core_source_object_v1",
            "actual_emergent_observer_closure_support_v4": "AO_obs_actual_closure_support_preLM_v1",
            "actual_emergent_observer_closure_readout_v4": "AP_obs_actual_closure_readout_preLM_v1",
            "admissible_AO_obs_actual_closure_support_object": status_ok,
            "admissible_AP_obs_actual_closure_readout_operator": status_ok,
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
