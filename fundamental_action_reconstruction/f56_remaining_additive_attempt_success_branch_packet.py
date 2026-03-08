#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f56_remaining_additive_attempt_success_branch_packet.json"
)
OUT_SUMMARY = (
    GENERATED / "f56_remaining_additive_attempt_success_branch_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n154 = load_json(
        "fundamental_action_reconstruction/generated/n154_next_constructive_move_reduced_to_one_explicit_additive_success_failure_branch_split_theorem_summary.json"
    )
    n155 = load_json(
        "fundamental_action_reconstruction/generated/n155_current_failure_branch_obstruction_theorem_for_additive_construction_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "n154_binary_branch_split_explicit",
            "actual": n154["theorem_result"]["binary_branch_split_explicit"],
            "expected": True,
            "meaning": "N154 already freezes the binary verdict split",
        },
        {
            "id": "n155_failure_branch_obstructed",
            "actual": n155["theorem_result"]["current_failure_branch_obstructed"],
            "expected": True,
            "meaning": "N155 already packages the current obstruction on failure_branch",
        },
        {
            "id": "n154_success_branch_name",
            "actual": n154["remaining_open_branches"][1],
            "expected": "future_success_branch_discharge_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "the success branch name is fixed",
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
        "stage": "F56",
        "lane": "remaining_success_branch_only",
        "goal": "freeze_the_success_branch_as_the_only_remaining_branch_to_attack_after_the_current_failure_branch_obstruction",
        "status": "F56_EXECUTED_REMAINING_ADDITIVE_ATTEMPT_SUCCESS_BRANCH_PACKET_NO_FALSE_PASS",
        "remaining_branch_to_attack": "explicit_success_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
        "basis": [
            "binary_split_already_frozen_by_N154",
            "failure_branch_currently_obstructed_by_N155",
            "no_third_verdict_branch_is_exported_at_this_scope",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F56",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "remaining_branch_to_attack": artifact["remaining_branch_to_attack"],
        "basis": artifact["basis"],
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
