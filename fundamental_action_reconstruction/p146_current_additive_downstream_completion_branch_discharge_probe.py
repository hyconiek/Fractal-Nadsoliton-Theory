#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "p146_current_additive_downstream_completion_branch_discharge_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f59 = load_json(
        "fundamental_action_reconstruction/generated/f59_last_additive_downstream_completion_branch_packet_summary.json"
    )
    n159 = load_json(
        "fundamental_action_reconstruction/generated/n159_current_additive_orientation_export_branch_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f59_last_remaining_lower_branch",
            "actual": f59["last_remaining_lower_branch"],
            "expected": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction_for_the_fixed_first_additive_attempt",
            "meaning": "F59 freezes the additive-specific downstream branch as the last remaining lower branch",
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
            "stage": "P146",
            "lane": "current_additive_downstream_completion_branch_discharge_only",
            "status": "P146_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADDITIVE_DOWNSTREAM_BRANCH_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P146",
            "lane": "current_additive_downstream_completion_branch_discharge_only",
            "goal": "test_whether_the_current_repo_already_exports_an_explicit_downstream_completion_branch_discharge_after_additive_specific_orientation_export_obstruction",
            "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_DOWNSTREAM_COMPLETION_BRANCH_DISCHARGE_FOR_THE_LAST_REMAINING_ADDITIVE_SPECIFIC_LOWER_BRANCH_AFTER_P146",
            "fixed_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "tested_downstream_branch": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction_for_the_fixed_first_additive_attempt",
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
