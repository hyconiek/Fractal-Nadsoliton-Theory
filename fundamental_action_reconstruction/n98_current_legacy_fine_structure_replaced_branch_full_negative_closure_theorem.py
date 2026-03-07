#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n98_current_legacy_fine_structure_replaced_branch_full_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n94 = load_json(
        "fundamental_action_reconstruction/generated/n94_current_strict_side_object_successor_branch_full_negative_closure_theorem_for_legacy_fine_structure_role_summary.json"
    )
    n97 = load_json(
        "fundamental_action_reconstruction/generated/n97_current_strict_side_method_successor_branch_full_negative_closure_theorem_for_legacy_fine_structure_role_summary.json"
    )

    checks_spec = [
        {
            "id": "n94_object_branch_closed",
            "actual": n94["theorem_result"]["object_successor_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N94 already closes the full object-successor branch negatively on the current repo state",
        },
        {
            "id": "n97_method_branch_closed",
            "actual": n97["theorem_result"]["method_successor_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N97 already closes the full method-successor branch negatively on the current repo state",
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
            "step": "N98",
            "status": "N98_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_FINE_STRUCTURE_REPLACED_BRANCH_STATE",
            "scope": "current_legacy_fine_structure_replaced_branch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N98",
            "status": "N98_DISCHARGED_CURRENT_LEGACY_FINE_STRUCTURE_REPLACED_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_fine_structure_replaced_branch_question_only",
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
