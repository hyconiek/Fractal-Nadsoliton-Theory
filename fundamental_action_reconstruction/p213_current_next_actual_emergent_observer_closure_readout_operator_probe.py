#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p213_current_next_actual_emergent_observer_closure_readout_operator_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n232 = load_json(
        GENERATED / "n232_current_next_actual_emergent_observer_closure_support_object_theorem_summary.json"
    )
    f125 = load_json(
        GENERATED / "f125_next_actual_emergent_observer_closure_readout_operator_packet_summary.json"
    )

    props = f125["observer_actual_closure_readout_properties"]
    response = f125["source_observer_actual_closure_readout"]
    checks_spec = [
        {
            "id": "actual_closure_support_is_admissible",
            "actual": n232["theorem_result"]["admissible_AO_obs_actual_closure_support_object"],
            "expected": True,
        },
        {
            "id": "derived_only_from_actual_closure_support_state",
            "actual": props["derived_only_from_actual_closure_support_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_actual_closure_readout_only",
            "actual": props["downstream_actual_closure_readout_only"]
            and not props["actual_emergent_observer_closure"],
            "expected": True,
        },
        {
            "id": "positive_commit_amplitude",
            "actual": response["positive_commit_amplitude"],
            "expected": True,
        },
        {
            "id": "zero_gap_channel",
            "actual": response["zero_gap_channel"],
            "expected": True,
        },
        {
            "id": "observer_information_deficit_downstream",
            "actual": props["observer_information_deficit_is_downstream_symptom"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe",
            "actual": props["kernel_split_safe"],
            "expected": True,
        },
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "stage": "P213",
            "lane": "current_next_actual_emergent_observer_closure_readout_only",
            "status": "P213_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P213",
            "lane": "current_next_actual_emergent_observer_closure_readout_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_FROM_AO_OBS_ACTUAL_CLOSURE_SUPPORT_PRELM_V1_AFTER_P213",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_actual_closure_support_object": "AO_obs_actual_closure_support_preLM_v1",
            "observer_actual_closure_readout_operator": "AP_obs_actual_closure_readout_preLM_v1",
            "checks": checks,
            "actual_closure_readout_exported": True,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
