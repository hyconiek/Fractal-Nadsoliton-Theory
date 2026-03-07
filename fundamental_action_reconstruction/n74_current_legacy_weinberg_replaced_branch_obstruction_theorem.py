#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n74_current_legacy_weinberg_replaced_branch_obstruction_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n73 = load_json(
        "fundamental_action_reconstruction/generated/n73_current_legacy_weinberg_retained_branch_full_negative_closure_theorem_summary.json"
    )
    p71 = load_json(
        "fundamental_action_reconstruction/generated/p71_strict_side_replaced_branch_probe_for_legacy_weinberg_role_summary.json"
    )

    checks_spec = [
        {
            "id": "n73_retained_branch_closed",
            "actual": n73["theorem_result"]["retained_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "the retained branch is already closed negatively on the current repo state",
        },
        {
            "id": "p71_probe_negative",
            "actual": p71["status"],
            "expected": "CURRENT_STRICT_SIDE_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_REPLACED_BRANCH_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P71",
            "meaning": "P71 confirms the replaced branch is still not exported on the current repo state",
        },
        {
            "id": "replaced_verdict_absent",
            "actual": p71["replaced_verdict_present"],
            "expected": False,
            "meaning": "explicit replaced-branch verdict is absent on the current repo state",
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
            "step": "N74",
            "status": "N74_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_WEINBERG_REPLACED_BRANCH_STATE",
            "scope": "current_legacy_weinberg_replaced_branch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N74",
            "status": "N74_DISCHARGED_CURRENT_LEGACY_WEINBERG_REPLACED_BRANCH_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_weinberg_replaced_branch_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "replaced_branch_present": False,
                "claim_specific_weinberg_frontier_still_open": True,
                "full_closure_pass": False,
            },
            "remaining_open_branch": [
                "explicit_strict_side_replaced_verdict_for_the_legacy_weinberg_angle_role_by_an_explicit_strict_successor_semantics"
            ],
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
