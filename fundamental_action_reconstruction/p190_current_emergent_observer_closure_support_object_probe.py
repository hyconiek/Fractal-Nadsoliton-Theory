#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p190_current_emergent_observer_closure_support_object_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n209 = load_json(
        "fundamental_action_reconstruction/generated/n209_current_emergent_observer_closure_fixed_point_test_theorem_summary.json"
    )
    f102 = load_json(
        "fundamental_action_reconstruction/generated/f102_first_emergent_observer_closure_support_object_packet_summary.json"
    )

    props = f102["observer_closure_support_properties"]
    response = f102["source_observer_closure_support"]
    checks_spec = [
        {
            "id": "observer_closure_fixed_point_test_is_admissible",
            "actual": n209["theorem_result"]["admissible_R_obs_closure_fixed_point_test"],
            "expected": True,
        },
        {
            "id": "derived_only_from_closure_fixed_point_state",
            "actual": props["derived_only_from_closure_fixed_point_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_closure_support_only",
            "actual": props["downstream_closure_support_only"]
            and not props["actual_emergent_observer_closure"],
            "expected": True,
        },
        {
            "id": "one_dimensional_closure_support_exported",
            "actual": response["one_dimensional_closure_support"],
            "expected": True,
        },
        {
            "id": "positive_closure_support_amplitude",
            "actual": response["positive_closure_support_amplitude"],
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
            "stage": "P190",
            "lane": "current_emergent_observer_closure_support_object_only",
            "status": "P190_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_EMERGENT_OBSERVER_CLOSURE_SUPPORT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P190",
            "lane": "current_emergent_observer_closure_support_object_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_FROM_R_OBS_CLOSURE_FIXED_POINT_TEST_PRELM_V1_AFTER_P190",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_closure_fixed_point_test": "R_obs_closure_fixed_point_test_preLM_v1",
            "observer_closure_support_map": "S_obs_closure_support_preLM_v1",
            "checks": checks,
            "closure_support_object_exported": True,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
