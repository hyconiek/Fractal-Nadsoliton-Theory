#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n161_current_fixed_first_additive_attempt_post_verdict_lower_branch_full_negative_closure_theorem_summary.json"
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
    n160 = load_json(
        "fundamental_action_reconstruction/generated/n160_current_additive_downstream_completion_branch_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n158_additive_admissibility_branch_closed",
            "actual": n158["theorem_result"][
                "current_additive_post_verdict_admissibility_branch_obstructed"
            ],
            "expected": True,
            "meaning": "N158 already closes the additive-specific admissibility branch negatively on the current repo state",
        },
        {
            "id": "n159_additive_orientation_branch_closed",
            "actual": n159["theorem_result"][
                "current_additive_orientation_export_branch_obstructed"
            ],
            "expected": True,
            "meaning": "N159 already closes the additive-specific orientation-export branch negatively on the current repo state",
        },
        {
            "id": "n160_additive_downstream_branch_closed",
            "actual": n160["theorem_result"][
                "current_additive_downstream_completion_branch_obstructed"
            ],
            "expected": True,
            "meaning": "N160 already closes the additive-specific downstream-completion branch negatively on the current repo state",
        },
        {
            "id": "n158_fixed_attempt_instance",
            "actual": n158["theorem_result"]["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N158 fixes the additive construction attempt under test",
        },
        {
            "id": "n159_fixed_attempt_instance",
            "actual": n159["theorem_result"]["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N159 fixes the same additive construction attempt under test",
        },
        {
            "id": "n160_fixed_attempt_instance",
            "actual": n160["theorem_result"]["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N160 fixes the same additive construction attempt under test",
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
            "step": "N161",
            "status": "N161_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIXED_FIRST_ADDITIVE_POST_VERDICT_LOWER_BRANCH_STATE",
            "scope": "current_fixed_first_additive_attempt_post_verdict_lower_branch_frontier_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N161",
            "status": "N161_DISCHARGED_CURRENT_FIXED_FIRST_ADDITIVE_ATTEMPT_POST_VERDICT_LOWER_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_fixed_first_additive_attempt_post_verdict_lower_branch_frontier_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "fixed_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "additive_admissibility_branch_closed_negatively_on_current_repo_state": True,
                "additive_orientation_export_branch_closed_negatively_on_current_repo_state": True,
                "additive_downstream_completion_branch_closed_negatively_on_current_repo_state": True,
                "fixed_first_additive_attempt_post_verdict_lower_branch_frontier_closed_negatively_on_current_repo_state": True,
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
