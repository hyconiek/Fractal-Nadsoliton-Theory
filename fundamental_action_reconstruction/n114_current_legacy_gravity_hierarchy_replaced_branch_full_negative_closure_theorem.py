#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n114_current_legacy_gravity_hierarchy_replaced_branch_full_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n110 = load_json(
        "fundamental_action_reconstruction/generated/n110_current_strict_side_object_successor_branch_full_negative_closure_theorem_for_legacy_gravity_hierarchy_role_summary.json"
    )
    n113 = load_json(
        "fundamental_action_reconstruction/generated/n113_current_strict_side_method_successor_branch_full_negative_closure_theorem_for_legacy_gravity_hierarchy_role_summary.json"
    )

    checks_spec = [
        {
            "id": "n110_object_branch_closed",
            "actual": n110["theorem_result"]["object_successor_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N110 already closes the full object-successor branch negatively on the current repo state",
        },
        {
            "id": "n113_method_branch_closed",
            "actual": n113["theorem_result"]["method_successor_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N113 already closes the full method-successor branch negatively on the current repo state",
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
            "step": "N114",
            "status": "N114_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_GRAVITY_HIERARCHY_REPLACED_BRANCH_STATE",
            "scope": "current_legacy_gravity_hierarchy_replaced_branch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N114",
            "status": "N114_DISCHARGED_CURRENT_LEGACY_GRAVITY_HIERARCHY_REPLACED_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_gravity_hierarchy_replaced_branch_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "object_successor_branch_closed_negatively_on_current_repo_state": True,
                "method_successor_branch_closed_negatively_on_current_repo_state": True,
                "replaced_branch_closed_negatively_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [],
            "hard_limits": [
                "no_proof_that_strict_side_successor_semantics_are_impossible_forever",
                "no_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT)


if __name__ == "__main__":
    main()
