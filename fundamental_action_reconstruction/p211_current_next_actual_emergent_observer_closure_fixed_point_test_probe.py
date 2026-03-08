#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p211_current_next_actual_emergent_observer_closure_fixed_point_test_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n230 = load_json(
        GENERATED / "n230_current_next_actual_emergent_observer_closure_realization_map_theorem_summary.json"
    )
    f123 = load_json(
        GENERATED / "f123_next_actual_emergent_observer_closure_fixed_point_test_packet_summary.json"
    )

    props = f123["observer_actual_closure_fixed_point_properties"]
    response = f123["source_observer_actual_closure_fixed_point"]
    checks_spec = [
        {
            "id": "actual_closure_realization_is_admissible",
            "actual": n230["theorem_result"]["admissible_AM_obs_actual_closure_realization"],
            "expected": True,
        },
        {
            "id": "derived_only_from_actual_closure_realization_map",
            "actual": props["derived_only_from_actual_closure_realization_map"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_actual_closure_fixed_point_only",
            "actual": props["downstream_actual_closure_fixed_point_only"]
            and not props["actual_emergent_observer_closure"],
            "expected": True,
        },
        {
            "id": "positive_actual_closure_fixed_point_amplitude",
            "actual": response["positive_actual_closure_fixed_point_amplitude"],
            "expected": True,
        },
        {
            "id": "one_dimensional_actual_closure_fixed_point",
            "actual": response["one_dimensional_actual_closure_fixed_point"],
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
            "stage": "P211",
            "lane": "current_next_actual_emergent_observer_closure_fixed_point_only",
            "status": "P211_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_FIXED_POINT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P211",
            "lane": "current_next_actual_emergent_observer_closure_fixed_point_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_FIXED_POINT_TEST_FROM_AM_OBS_ACTUAL_CLOSURE_REALIZATION_PRELM_V1_AFTER_P211",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_actual_closure_realization_operator": "AM_obs_actual_closure_realization_preLM_v1",
            "observer_actual_closure_fixed_point_test": "AN_obs_actual_closure_fixed_point_test_preLM_v1",
            "checks": checks,
            "actual_closure_fixed_point_exported": True,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
