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
    / "p142_current_failure_verdict_discharge_probe_for_additive_construction_attempt.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p142_current_failure_verdict_discharge_probe_for_additive_construction_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f55 = load_json(
        "fundamental_action_reconstruction/generated/f55_first_conservative_additive_attempt_failure_branch_packet_summary.json"
    )
    n154 = load_json(
        "fundamental_action_reconstruction/generated/n154_next_constructive_move_reduced_to_one_explicit_additive_success_failure_branch_split_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f55_first_branch_to_attack",
            "actual": f55["first_branch_to_attack"],
            "expected": "explicit_failure_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "F55 freezes the additive failure branch as first branch to attack",
        },
        {
            "id": "n154_binary_branch_split_explicit",
            "actual": n154["theorem_result"]["binary_branch_split_explicit"],
            "expected": True,
            "meaning": "N154 already fixes the binary verdict split",
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
        "stage": "P142",
        "lane": "current_failure_verdict_discharge_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_failure_verdict_for_the_fixed_first_future_additive_construction_attempt",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_FAILURE_VERDICT_DISCHARGE_FOR_THE_FIXED_FIRST_ADDITIVE_CONSTRUCTION_ATTEMPT_AFTER_P142",
        "fixed_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
        "tested_failure_branch": "explicit_failure_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P142",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "fixed_attempt_instance": artifact["fixed_attempt_instance"],
        "tested_failure_branch": artifact["tested_failure_branch"],
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
