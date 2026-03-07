#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n83_current_legacy_weinberg_full_claim_specific_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n73 = load_json(
        "fundamental_action_reconstruction/generated/n73_current_legacy_weinberg_retained_branch_full_negative_closure_theorem_summary.json"
    )
    n82 = load_json(
        "fundamental_action_reconstruction/generated/n82_current_legacy_weinberg_replaced_branch_full_negative_closure_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n73_retained_branch_closed",
            "actual": n73["theorem_result"]["retained_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N73 already closes the retained branch negatively on the current repo state",
        },
        {
            "id": "n82_replaced_branch_closed",
            "actual": n82["theorem_result"]["replaced_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N82 already closes the replaced branch negatively on the current repo state",
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
            "step": "N83",
            "status": "N83_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_WEINBERG_FULL_CLAIM_SPECIFIC_STATE",
            "scope": "current_legacy_weinberg_full_claim_specific_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N83",
            "status": "N83_DISCHARGED_CURRENT_LEGACY_WEINBERG_FULL_CLAIM_SPECIFIC_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_weinberg_full_claim_specific_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "retained_branch_closed_negatively_on_current_repo_state": True,
                "replaced_branch_closed_negatively_on_current_repo_state": True,
                "legacy_weinberg_claim_specific_frontier_closed_negatively_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [],
            "hard_limits": [
                "no_proof_that_legacy_weinberg_transfer_is_impossible_forever",
                "no_proof_that_future_strict_side_successor_semantics_cannot_exist",
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
