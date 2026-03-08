#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f61_current_selector_construction_stop_condition_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
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
            "id": "n149_constructive_exhaustion",
            "actual": n149["theorem_result"][
                "constructive_selector_frontier_exhausted_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N149 already exhausts the constructive selector frontier on current exports",
        },
        {
            "id": "n162_fixed_additive_attempt_negative",
            "actual": n162["theorem_result"][
                "fixed_first_additive_construction_attempt_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N162 already closes the fixed first additive attempt negatively",
        },
        {
            "id": "n163_observer_deficit_downstream",
            "actual": n163["theorem_result"][
                "observer_information_deficit_downstream_symptom_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N163 already classifies observer information deficit as downstream symptom",
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
            "stage": "F61",
            "lane": "current_selector_construction_stop_condition_only",
            "status": "F61_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STOP_CONDITION_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F61",
            "lane": "current_selector_construction_stop_condition_only",
            "goal": "freeze_the_current_repo_state_stop_condition_for_the_current_selector_construction_lane",
            "status": "F61_EXECUTED_CURRENT_SELECTOR_CONSTRUCTION_STOP_CONDITION_PACKET_NO_FALSE_PASS",
            "checks": checks,
            "stop_condition": {
                "current_selector_construction_lane_stop_condition_active": True,
                "remaining_admissible_move_class": "future_genuinely_additive_upstream_source_work_only",
            },
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
