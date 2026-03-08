#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n162_current_fixed_first_additive_construction_attempt_full_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n157 = load_json(
        "fundamental_action_reconstruction/generated/n157_current_additive_construction_attempt_verdict_layer_full_negative_closure_theorem_summary.json"
    )
    n161 = load_json(
        "fundamental_action_reconstruction/generated/n161_current_fixed_first_additive_attempt_post_verdict_lower_branch_full_negative_closure_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n157_verdict_layer_closed",
            "actual": n157["theorem_result"][
                "additive_construction_attempt_verdict_layer_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N157 already closes the fixed first additive attempt verdict layer negatively on the current repo state",
        },
        {
            "id": "n161_lower_branch_frontier_closed",
            "actual": n161["theorem_result"][
                "fixed_first_additive_attempt_post_verdict_lower_branch_frontier_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N161 already closes the fixed first additive attempt post-verdict lower-branch frontier negatively on the current repo state",
        },
        {
            "id": "n157_fixed_attempt_instance",
            "actual": n157["theorem_result"]["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N157 fixes the additive construction attempt under test",
        },
        {
            "id": "n161_fixed_attempt_instance",
            "actual": n161["theorem_result"]["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N161 fixes the same additive construction attempt under test",
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
            "step": "N162",
            "status": "N162_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIXED_FIRST_ADDITIVE_ATTEMPT_STATE",
            "scope": "current_fixed_first_additive_construction_attempt_as_a_whole",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N162",
            "status": "N162_DISCHARGED_CURRENT_FIXED_FIRST_ADDITIVE_CONSTRUCTION_ATTEMPT_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_fixed_first_additive_construction_attempt_as_a_whole",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "fixed_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "verdict_layer_closed_negatively_on_current_repo_state": True,
                "post_verdict_lower_branch_frontier_closed_negatively_on_current_repo_state": True,
                "fixed_first_additive_construction_attempt_closed_negatively_on_current_repo_state": True,
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
