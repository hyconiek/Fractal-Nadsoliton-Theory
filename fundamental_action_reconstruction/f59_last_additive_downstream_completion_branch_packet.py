#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "f59_last_additive_downstream_completion_branch_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n158 = load_json(
        "fundamental_action_reconstruction/generated/n158_current_additive_post_verdict_admissibility_branch_obstruction_theorem_summary.json"
    )
    n159 = load_json(
        "fundamental_action_reconstruction/generated/n159_current_additive_orientation_export_branch_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n158_additive_admissibility_branch_closed",
            "actual": n158["theorem_result"]["current_additive_post_verdict_admissibility_branch_obstructed"],
            "expected": True,
            "meaning": "N158 already closes the additive-specific admissibility branch negatively on the current repo state",
        },
        {
            "id": "n159_additive_orientation_branch_closed",
            "actual": n159["theorem_result"]["current_additive_orientation_export_branch_obstructed"],
            "expected": True,
            "meaning": "N159 already closes the additive-specific orientation-export branch negatively on the current repo state",
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
            "stage": "F59",
            "lane": "last_additive_downstream_branch_only",
            "status": "F59_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADDITIVE_DOWNSTREAM_BRANCH_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F59",
            "lane": "last_additive_downstream_branch_only",
            "goal": "freeze_the_additive_specific_downstream_completion_branch_as_the_only_remaining_lower_branch",
            "status": "F59_EXECUTED_LAST_ADDITIVE_DOWNSTREAM_COMPLETION_BRANCH_PACKET_NO_FALSE_PASS",
            "last_remaining_lower_branch": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction_for_the_fixed_first_additive_attempt",
            "basis": [
                "additive_specific_admissibility_branch_already_closed_by_N158",
                "additive_specific_orientation_export_branch_already_closed_by_N159",
                "downstream_B_sel_R_sel_O_sel_still_comes_after_source_and_orientation_stages",
                "preferred_order_remains_nadsoliton_to_light_to_matter_to_emergent_observer",
                "observer_information_deficit_remains_downstream_of_source_and_export_stages",
            ],
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT)


if __name__ == "__main__":
    main()
