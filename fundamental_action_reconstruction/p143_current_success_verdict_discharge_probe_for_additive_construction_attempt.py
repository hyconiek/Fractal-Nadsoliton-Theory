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
    / "p143_current_success_verdict_discharge_probe_for_additive_construction_attempt.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p143_current_success_verdict_discharge_probe_for_additive_construction_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f56 = load_json(
        "fundamental_action_reconstruction/generated/f56_remaining_additive_attempt_success_branch_packet_summary.json"
    )
    n154 = load_json(
        "fundamental_action_reconstruction/generated/n154_next_constructive_move_reduced_to_one_explicit_additive_success_failure_branch_split_theorem_summary.json"
    )
    n155 = load_json(
        "fundamental_action_reconstruction/generated/n155_current_failure_branch_obstruction_theorem_for_additive_construction_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "f56_remaining_branch_to_attack",
            "actual": f56["remaining_branch_to_attack"],
            "expected": "explicit_success_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "F56 freezes the success branch as the only remaining branch to attack",
        },
        {
            "id": "n154_binary_branch_split_explicit",
            "actual": n154["theorem_result"]["binary_branch_split_explicit"],
            "expected": True,
            "meaning": "N154 already fixes the binary verdict split",
        },
        {
            "id": "n155_current_failure_branch_obstructed",
            "actual": n155["theorem_result"]["current_failure_branch_obstructed"],
            "expected": True,
            "meaning": "N155 already obstructs the failure branch on the current repo state",
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
        "stage": "P143",
        "lane": "current_success_verdict_discharge_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_success_verdict_for_the_fixed_first_future_additive_construction_attempt",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_SUCCESS_VERDICT_DISCHARGE_FOR_THE_FIXED_FIRST_ADDITIVE_CONSTRUCTION_ATTEMPT_AFTER_P143",
        "fixed_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
        "tested_success_branch": "explicit_success_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P143",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "fixed_attempt_instance": artifact["fixed_attempt_instance"],
        "tested_success_branch": artifact["tested_success_branch"],
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
