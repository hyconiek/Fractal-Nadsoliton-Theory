#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p189_current_emergent_observer_closure_fixed_point_test_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n208 = load_json(
        "fundamental_action_reconstruction/generated/n208_current_emergent_observer_closure_realization_map_theorem_summary.json"
    )
    f101 = load_json(
        "fundamental_action_reconstruction/generated/f101_first_emergent_observer_closure_fixed_point_test_packet_summary.json"
    )

    props = f101["observer_closure_fixed_point_properties"]
    response = f101["source_observer_closure_fixed_point"]
    checks_spec = [
        {
            "id": "observer_closure_realization_is_admissible",
            "actual": n208["theorem_result"]["admissible_Q_obs_closure_realization"],
            "expected": True,
        },
        {
            "id": "derived_only_from_closure_realization_state",
            "actual": props["derived_only_from_closure_realization_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_closure_fixed_point_only",
            "actual": props["downstream_closure_fixed_point_only"]
            and not props["actual_emergent_observer_closure"],
            "expected": True,
        },
        {
            "id": "one_dimensional_closure_fixed_point_exported",
            "actual": response["one_dimensional_closure_fixed_point"],
            "expected": True,
        },
        {
            "id": "positive_closure_fixed_point_amplitude",
            "actual": response["positive_closure_fixed_point_amplitude"],
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
            "stage": "P189",
            "lane": "current_emergent_observer_closure_fixed_point_test_only",
            "status": "P189_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_EMERGENT_OBSERVER_CLOSURE_FIXED_POINT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P189",
            "lane": "current_emergent_observer_closure_fixed_point_test_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_CLOSURE_FIXED_POINT_TEST_FROM_Q_OBS_CLOSURE_REALIZATION_PRELM_V1_AFTER_P189",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_closure_realization_map": "Q_obs_closure_realization_preLM_v1",
            "observer_closure_fixed_point_test": "R_obs_closure_fixed_point_test_preLM_v1",
            "checks": checks,
            "closure_fixed_point_test_exported": True,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
