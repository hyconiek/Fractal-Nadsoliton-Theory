#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p201_current_actual_emergent_observer_closure_readout_operator_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n220 = load_json(GENERATED / "n220_current_actual_emergent_observer_closure_support_object_theorem_summary.json")
    f113 = load_json(GENERATED / "f113_first_actual_emergent_observer_closure_readout_operator_packet_summary.json")

    checks = [
        {
            "id": "observer_actual_closure_support_is_admissible",
            "actual": n220["theorem_result"]["admissible_AC_obs_actual_closure_support_object"],
            "expected": True,
        },
        {
            "id": "derived_only_from_actual_closure_support",
            "actual": f113["observer_actual_closure_readout_properties"]["derived_only_from_actual_closure_support"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": f113["observer_actual_closure_readout_properties"]["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_actual_closure_readout_only",
            "actual": f113["observer_actual_closure_readout_properties"]["downstream_actual_closure_readout_only"],
            "expected": True,
        },
        {
            "id": "positive_actual_closure_commit",
            "actual": f113["source_observer_actual_closure_readout"]["positive_actual_closure_commit"],
            "expected": True,
        },
        {
            "id": "zero_actual_closure_gap",
            "actual": f113["source_observer_actual_closure_readout"]["zero_actual_closure_gap"],
            "expected": True,
        },
        {
            "id": "observer_information_deficit_downstream",
            "actual": f113["observer_actual_closure_readout_properties"]["observer_information_deficit_is_downstream_symptom"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe",
            "actual": f113["observer_actual_closure_readout_properties"]["kernel_split_safe"],
            "expected": True,
        },
    ]

    for check in checks:
        check["pass"] = check["actual"] == check["expected"]

    status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_FROM_AC_OBS_ACTUAL_CLOSURE_SUPPORT_PRELM_V1_AFTER_P201"
        if all(check["pass"] for check in checks)
        else "CURRENT_REPO_DOES_NOT_EXPORT_AN_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_FROM_AC_OBS_ACTUAL_CLOSURE_SUPPORT_PRELM_V1_AFTER_P201"
    )

    summary = {
        "stage": "P201",
        "lane": "current_actual_emergent_observer_closure_readout_only",
        "status": status,
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_support_object": "AC_obs_actual_closure_support_preLM_v1",
        "observer_actual_closure_readout_operator": "AD_obs_actual_closure_readout_preLM_v1",
        "checks": checks,
        "actual_closure_readout_exported": all(check["pass"] for check in checks),
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
