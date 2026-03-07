#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n81_current_strict_side_method_successor_branch_full_negative_closure_theorem_for_legacy_weinberg_role_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n80 = load_json(
        "fundamental_action_reconstruction/generated/n80_current_strict_side_textual_method_successor_obstruction_theorem_for_legacy_weinberg_role_summary.json"
    )
    p78 = load_json(
        "fundamental_action_reconstruction/generated/p78_strict_side_method_lineage_upgrade_probe_for_legacy_weinberg_role_summary.json"
    )

    checks_spec = [
        {
            "id": "n80_textual_method_successor_closed",
            "actual": n80["theorem_result"]["textual_method_successor_subbranch_present"],
            "expected": False,
            "meaning": "N80 already keeps the textual method-successor sub-branch absent",
        },
        {
            "id": "p78_method_lineage_upgrade_absent",
            "actual": p78["method_lineage_upgrade_verdict_present"],
            "expected": False,
            "meaning": "P78 keeps the method-lineage-upgrade sub-branch absent",
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
            "step": "N81",
            "status": "N81_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_SIDE_METHOD_SUCCESSOR_BRANCH_STATE",
            "scope": "current_strict_side_method_successor_branch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N81",
            "status": "N81_DISCHARGED_CURRENT_STRICT_SIDE_METHOD_SUCCESSOR_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM_FOR_LEGACY_WEINBERG_ROLE_NO_FALSE_PASS",
            "scope": "current_strict_side_method_successor_branch_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "textual_method_successor_subbranch_present": False,
                "method_lineage_upgrade_subbranch_present": False,
                "method_successor_branch_closed_negatively_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [],
            "hard_limits": [
                "no_proof_that_the_replaced_branch_is_impossible_forever",
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
