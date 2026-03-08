#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "p145_current_additive_orientation_export_branch_discharge_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f58 = load_json(
        "fundamental_action_reconstruction/generated/f58_first_post_additive_admissibility_orientation_export_branch_packet_summary.json"
    )
    n158 = load_json(
        "fundamental_action_reconstruction/generated/n158_current_additive_post_verdict_admissibility_branch_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f58_first_remaining_lower_branch_to_attack",
            "actual": f58["first_remaining_lower_branch_to_attack"],
            "expected": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object_for_the_fixed_first_additive_attempt",
            "meaning": "F58 freezes the additive-specific E_orient branch as the first remaining lower branch",
        },
        {
            "id": "n158_additive_admissibility_branch_closed",
            "actual": n158["theorem_result"]["current_additive_post_verdict_admissibility_branch_obstructed"],
            "expected": True,
            "meaning": "N158 already closes the additive-specific admissibility branch negatively on the current repo state",
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
            "stage": "P145",
            "lane": "current_additive_orientation_export_branch_discharge_only",
            "status": "P145_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADDITIVE_ORIENTATION_BRANCH_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P145",
            "lane": "current_additive_orientation_export_branch_discharge_only",
            "goal": "test_whether_the_current_repo_already_exports_an_explicit_orientation_export_branch_discharge_after_additive_admissibility_obstruction",
            "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_ORIENTATION_EXPORT_BRANCH_DISCHARGE_FOR_THE_FIRST_REMAINING_FUTURE_ADDITIVE_E_ORIENT_BRANCH_AFTER_P145",
            "fixed_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "tested_orientation_branch": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object_for_the_fixed_first_additive_attempt",
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
