#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f62_current_selector_construction_handoff_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n149 = load_json(
        "fundamental_action_reconstruction/generated/n149_current_repo_constructive_selector_frontier_exhaustion_theorem_summary.json"
    )
    n163 = load_json(
        "fundamental_action_reconstruction/generated/n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
    )
    n164 = load_json(
        "fundamental_action_reconstruction/generated/n164_current_selector_construction_stop_condition_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n149_exhausted",
            "actual": n149["theorem_result"][
                "constructive_selector_frontier_exhausted_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N149 already exhausts the constructive selector frontier",
        },
        {
            "id": "n163_downstream_observer_deficit",
            "actual": n163["theorem_result"][
                "observer_information_deficit_downstream_symptom_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N163 already keeps observer deficit downstream",
        },
        {
            "id": "n164_lane_stopped",
            "actual": n164["theorem_result"][
                "current_selector_construction_lane_should_be_treated_as_stopped"
            ],
            "expected": True,
            "meaning": "N164 already stops the present selector-construction lane",
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
            "stage": "F62",
            "lane": "current_selector_construction_handoff_only",
            "status": "F62_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_HANDOFF_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F62",
            "lane": "current_selector_construction_handoff_only",
            "goal": "freeze_the_current_repo_state_handoff_from_the_stopped_selector_construction_lane_to_genuinely_additive_future_upstream_source_work",
            "status": "F62_EXECUTED_CURRENT_SELECTOR_CONSTRUCTION_HANDOFF_PACKET_NO_FALSE_PASS",
            "checks": checks,
            "handoff": {
                "current_selector_construction_handoff_active": True,
                "handoff_target": "future_genuinely_additive_upstream_source_work_only",
            },
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
