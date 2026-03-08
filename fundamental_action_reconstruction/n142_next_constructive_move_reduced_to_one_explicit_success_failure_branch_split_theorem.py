#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n142_next_constructive_move_reduced_to_one_explicit_success_failure_branch_split_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p131 = load_json(
        "fundamental_action_reconstruction/generated/p131_first_future_constructed_source_object_realization_verdict_branch_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p131_binary_branch_split_reduced",
            "actual": p131["target_state"][
                "next_constructive_move_reduced_to_binary_branch_split"
            ],
            "expected": True,
            "meaning": "P131 already proves the next constructive move is reduced to one explicit binary branch split",
        },
        {
            "id": "success_branch_name",
            "actual": p131["target_state"]["verdict_branches"]["success_branch"],
            "expected": "explicit_success_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0",
            "meaning": "the success branch is explicitly named",
        },
        {
            "id": "failure_branch_name",
            "actual": p131["target_state"]["verdict_branches"]["failure_branch"],
            "expected": "explicit_failure_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0",
            "meaning": "the failure branch is explicitly named",
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
            "step": "N142",
            "status": "N142_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SUCCESS_FAILURE_BRANCH_SPLIT_STATE",
            "scope": "first_future_constructed_source_object_realization_verdict_branch_refinement_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N142",
            "status": "N142_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_EXPLICIT_SUCCESS_FAILURE_BRANCH_SPLIT_THEOREM_NO_FALSE_PASS",
            "scope": "first_future_constructed_source_object_realization_verdict_branch_refinement_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_constructive_move_reduced_to_binary_branch_split": True,
                "verdict_target": p131["target_state"]["verdict_target"],
                "verdict_branches": p131["target_state"]["verdict_branches"],
                "later_open_branches": p131["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_success_verdict_discharge_for_S_sel_int_new_object_constructed_realization_attempt_v0",
                "future_failure_verdict_discharge_for_S_sel_int_new_object_constructed_realization_attempt_v0",
                "future_admissibility_test_of_a_future_constructed_source_object",
                "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
                "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            ],
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
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
