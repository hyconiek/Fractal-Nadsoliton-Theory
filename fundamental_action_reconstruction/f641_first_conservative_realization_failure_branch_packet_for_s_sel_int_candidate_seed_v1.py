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
    / "f641_first_conservative_realization_failure_branch_packet_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f641_first_conservative_realization_failure_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n531 = load_json(
        "fundamental_action_reconstruction/generated/n531_next_constructive_move_reduced_to_one_explicit_success_failure_branch_split_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    checks_spec = [
        {
            "id": "n531_binary_branch_split_reduced",
            "actual": n531["theorem_result"][
                "next_constructive_move_reduced_to_binary_branch_split"
            ],
            "expected": True,
            "meaning": "N531 already freezes the explicit binary success/failure verdict split (v1)",
        },
        {
            "id": "success_branch_name",
            "actual": n531["theorem_result"]["verdict_branches"]["success_branch"],
            "expected": "explicit_success_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v1",
            "meaning": "the success branch name is fixed (v1)",
        },
        {
            "id": "failure_branch_name",
            "actual": n531["theorem_result"]["verdict_branches"]["failure_branch"],
            "expected": "explicit_failure_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v1",
            "meaning": "the failure branch name is fixed (v1)",
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
        "stage": "F641",
        "lane": "conservative_failure_branch_first_only",
        "goal": "freeze_the_failure_branch_as_the_first_branch_to_attack_under_no_false_pass_ordering_for_seed_v1",
        "status": "F641_EXECUTED_FIRST_CONSERVATIVE_REALIZATION_FAILURE_BRANCH_PACKET_FOR_SEED_V1_NO_FALSE_PASS",
        "first_branch_to_attack": "explicit_failure_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v1",
        "branch_ordering_basis": [
            "binary_split_already_frozen_by_N531",
            "no_realized_constructed_source_object",
            "no_admissible_S_sel_int",
            "failure_first_is_methodologically_more_conservative_than_success_first",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F641",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "first_branch_to_attack": artifact["first_branch_to_attack"],
        "branch_ordering_basis": artifact["branch_ordering_basis"],
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

