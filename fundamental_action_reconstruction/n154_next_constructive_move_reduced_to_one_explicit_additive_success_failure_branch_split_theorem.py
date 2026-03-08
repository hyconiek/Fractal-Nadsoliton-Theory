#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n154_next_constructive_move_reduced_to_one_explicit_additive_success_failure_branch_split_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p141 = load_json(
        "fundamental_action_reconstruction/generated/p141_first_future_additive_construction_attempt_verdict_branch_probe_summary.json"
    )
    n153 = load_json(
        "fundamental_action_reconstruction/generated/n153_next_constructive_move_reduced_to_one_first_additive_construction_attempt_verdict_target_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n153_verdict_target",
            "actual": n153["theorem_result"]["explicit_future_additive_attempt_verdict_target"],
            "expected": "success_or_failure_verdict(construct_attempt_v1(S_sel_int_additive_attempt_target_v1))",
            "meaning": "N153 freezes one fixed additive-construction verdict target",
        },
        {
            "id": "p141_remaining_open_branches",
            "actual": p141["target_state"]["remaining_open_branches"],
            "expected": [
                "future_failure_branch_discharge_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "future_success_branch_discharge_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)"
            ],
            "meaning": "P141 reduces the next move to one explicit binary branch split",
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
            "step": "N154",
            "status": "N154_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADDITIVE_SUCCESS_FAILURE_BRANCH_SPLIT_STATE",
            "scope": "next_constructive_move_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N154",
            "status": "N154_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_EXPLICIT_ADDITIVE_SUCCESS_FAILURE_BRANCH_SPLIT_THEOREM_NO_FALSE_PASS",
            "scope": "next_constructive_move_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "fixed_verdict_target": "success_or_failure_verdict(construct_attempt_v1(S_sel_int_additive_attempt_target_v1))",
                "binary_branch_split_explicit": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_failure_branch_discharge_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "future_success_branch_discharge_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)"
            ],
            "hard_limits": [
                "success_verdict_not_yet_exported",
                "failure_verdict_not_yet_exported",
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
