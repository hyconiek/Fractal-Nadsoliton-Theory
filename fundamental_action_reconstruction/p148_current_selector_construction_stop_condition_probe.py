#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p148_current_selector_construction_stop_condition_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f61 = load_json(
        "fundamental_action_reconstruction/generated/f61_current_selector_construction_stop_condition_packet_summary.json"
    )
    n149 = load_json(
        "fundamental_action_reconstruction/generated/n149_current_repo_constructive_selector_frontier_exhaustion_theorem_summary.json"
    )
    n162 = load_json(
        "fundamental_action_reconstruction/generated/n162_current_fixed_first_additive_construction_attempt_full_negative_closure_theorem_summary.json"
    )
    n163 = load_json(
        "fundamental_action_reconstruction/generated/n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f61_stop_condition_active",
            "actual": f61["stop_condition"]["current_selector_construction_lane_stop_condition_active"],
            "expected": True,
            "meaning": "F61 already freezes the stop condition as active",
        },
        {
            "id": "n149_exhausted",
            "actual": n149["theorem_result"][
                "constructive_selector_frontier_exhausted_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N149 already exhausts current constructive selector frontier",
        },
        {
            "id": "n162_fixed_additive_negative",
            "actual": n162["theorem_result"][
                "fixed_first_additive_construction_attempt_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N162 already closes the fixed first additive attempt negatively",
        },
        {
            "id": "n163_downstream_observer_deficit",
            "actual": n163["theorem_result"][
                "observer_information_deficit_downstream_symptom_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N163 already keeps observer deficit downstream",
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
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "stage": "P148",
            "lane": "current_selector_construction_stop_condition_probe_only",
            "status": "P148_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STOP_CONDITION_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P148",
            "lane": "current_selector_construction_stop_condition_probe_only",
            "goal": "test_whether_the_current_repo_supports_stopping_the_present_selector_construction_lane_and_leaving_only_genuinely_additive_future_upstream_source_work",
            "status": "CURRENT_REPO_SUPPORTS_THE_CONCLUSION_THAT_THE_PRESENT_SELECTOR_CONSTRUCTION_LANE_SHOULD_BE_TREATED_AS_STOPPED_AND_ONLY_GENUINELY_ADDITIVE_FUTURE_UPSTREAM_SOURCE_WORK_REMAINS_AFTER_P148",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
