#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n531_next_constructive_move_reduced_to_one_explicit_success_failure_branch_split_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p640 = load_json(
        "fundamental_action_reconstruction/generated/p640_first_future_constructed_source_object_realization_verdict_branch_probe_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    checks_spec = [
        {
            "id": "p640_binary_branch_split_reduced",
            "actual": p640["target_state"][
                "next_constructive_move_reduced_to_binary_branch_split"
            ],
            "expected": True,
            "meaning": "P640 proves the next constructive move is reduced to one explicit binary branch split (v1)",
        },
        {
            "id": "success_branch_name",
            "actual": p640["target_state"]["verdict_branches"]["success_branch"],
            "expected": "explicit_success_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v1",
            "meaning": "the success branch is explicitly named (v1)",
        },
        {
            "id": "failure_branch_name",
            "actual": p640["target_state"]["verdict_branches"]["failure_branch"],
            "expected": "explicit_failure_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v1",
            "meaning": "the failure branch is explicitly named (v1)",
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
            "step": "N531",
            "status": "N531_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SUCCESS_FAILURE_BRANCH_SPLIT_STATE_FOR_SEED_V1",
            "scope": "first_future_constructed_source_object_realization_verdict_branch_refinement_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
    else:
        summary = {
            "step": "N531",
            "status": "N531_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_EXPLICIT_SUCCESS_FAILURE_BRANCH_SPLIT_THEOREM_FOR_SEED_V1_NO_FALSE_PASS",
            "scope": "first_future_constructed_source_object_realization_verdict_branch_refinement_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_constructive_move_reduced_to_binary_branch_split": True,
                "verdict_target": p640["target_state"]["verdict_target"],
                "verdict_branches": p640["target_state"]["verdict_branches"],
                "later_open_branches": p640["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": p640["target_state"]["later_open_branches"],
            "hard_limits": [
                "success_branch_not_yet_discharged",
                "failure_branch_not_yet_discharged",
                "constructed_source_object_not_yet_exported",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

