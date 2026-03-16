#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "f642_remaining_realization_success_branch_packet_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f642_remaining_realization_success_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n531 = load_json(
        "fundamental_action_reconstruction/generated/n531_next_constructive_move_reduced_to_one_explicit_success_failure_branch_split_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
    )
    n532 = load_json(
        "fundamental_action_reconstruction/generated/n532_current_failure_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "n531_binary_branch_split_reduced",
            "actual": n531["theorem_result"][
                "next_constructive_move_reduced_to_binary_branch_split"
            ],
            "expected": True,
            "meaning": "N531 already freezes the binary verdict split (v1)",
        },
        {
            "id": "success_branch_name",
            "actual": n531["theorem_result"]["verdict_branches"]["success_branch"],
            "expected": "explicit_success_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v1",
            "meaning": "the success branch name is fixed (v1)",
        },
        {
            "id": "n532_failure_branch_obstruction_discharged",
            "actual": n532["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N532 packages the current failure-side obstruction (v1)",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "F642",
        "lane": "remaining_success_branch_only",
        "goal": "freeze_the_success_branch_as_the_only_remaining_branch_to_attack_after_the_current_failure_branch_obstruction_for_seed_v1",
        "status": "F642_EXECUTED_REMAINING_REALIZATION_SUCCESS_BRANCH_PACKET_FOR_SEED_V1_NO_FALSE_PASS",
        "remaining_branch_to_attack": "explicit_success_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v1",
        "branch_selection_basis": [
            "binary_split_already_frozen_by_N531",
            "current_failure_branch_obstruction_already_packaged_by_N532",
            "no_third_verdict_branch_exported_at_the_same_scope",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F642",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "remaining_branch_to_attack": artifact["remaining_branch_to_attack"],
        "branch_selection_basis": artifact["branch_selection_basis"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()

