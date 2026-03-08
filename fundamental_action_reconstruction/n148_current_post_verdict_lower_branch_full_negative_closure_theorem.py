#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n148_current_post_verdict_lower_branch_full_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n145 = load_json(
        "fundamental_action_reconstruction/generated/n145_current_admissibility_branch_obstruction_theorem_for_future_constructed_source_object_summary.json"
    )
    n146 = load_json(
        "fundamental_action_reconstruction/generated/n146_current_orientation_export_branch_obstruction_theorem_summary.json"
    )
    n147 = load_json(
        "fundamental_action_reconstruction/generated/n147_current_downstream_completion_branch_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n145_admissibility_branch_closed",
            "actual": n145["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N145 already closes the admissibility branch negatively on the current repo state",
        },
        {
            "id": "n146_orientation_export_branch_closed",
            "actual": n146["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N146 already closes the orientation-export branch negatively on the current repo state",
        },
        {
            "id": "n147_downstream_branch_closed",
            "actual": n147["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N147 already closes the downstream-completion branch negatively on the current repo state",
        },
        {
            "id": "n147_remaining_open_branches_empty",
            "actual": n147["remaining_open_branches"],
            "expected": [],
            "meaning": "after N147 the lower-branch list is exhausted on the current repo state",
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
            "step": "N148",
            "status": "N148_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_POST_VERDICT_LOWER_BRANCH_STATE",
            "scope": "current_post_verdict_lower_branch_frontier_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N148",
            "status": "N148_DISCHARGED_CURRENT_POST_VERDICT_LOWER_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_post_verdict_lower_branch_frontier_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "admissibility_branch_closed_negatively_on_current_repo_state": True,
                "orientation_export_branch_closed_negatively_on_current_repo_state": True,
                "downstream_completion_branch_closed_negatively_on_current_repo_state": True,
                "post_verdict_lower_branch_frontier_closed_negatively_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [],
            "hard_limits": [
                "no_proof_that_future_additive_source_object_construction_is_impossible_forever",
                "constructed_source_object_not_yet_exported",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
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
