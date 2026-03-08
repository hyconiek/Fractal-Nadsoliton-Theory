#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n153_next_constructive_move_reduced_to_one_first_additive_construction_attempt_verdict_target_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p140 = load_json(
        "fundamental_action_reconstruction/generated/p140_first_future_additive_construction_attempt_verdict_target_probe_summary.json"
    )
    n152 = load_json(
        "fundamental_action_reconstruction/generated/n152_only_remaining_positive_move_class_reduced_to_one_first_additive_construction_attempt_instance_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n152_attempt_instance",
            "actual": n152["theorem_result"]["explicit_future_additive_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N152 fixes the future additive construction-attempt instance",
        },
        {
            "id": "p140_verdict_target",
            "actual": p140["target_state"]["explicit_future_additive_attempt_verdict_target"],
            "expected": "success_or_failure_verdict(construct_attempt_v1(S_sel_int_additive_attempt_target_v1))",
            "meaning": "P140 reduces the next constructive move to one verdict target",
        },
        {
            "id": "p140_single_remaining_branch",
            "actual": p140["target_state"]["remaining_open_branches"],
            "expected": [
                "future_explicit_success_failure_branch_split_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)"
            ],
            "meaning": "P140 leaves only one future branch-split move",
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
            "step": "N153",
            "status": "N153_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADDITIVE_ATTEMPT_VERDICT_TARGET_REDUCTION_STATE",
            "scope": "next_constructive_move_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N153",
            "status": "N153_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_FIRST_ADDITIVE_CONSTRUCTION_ATTEMPT_VERDICT_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "next_constructive_move_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "explicit_future_additive_attempt_verdict_target": "success_or_failure_verdict(construct_attempt_v1(S_sel_int_additive_attempt_target_v1))",
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_explicit_success_failure_branch_split_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)"
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
