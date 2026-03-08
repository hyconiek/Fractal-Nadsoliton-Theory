#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f55_first_conservative_additive_attempt_failure_branch_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f55_first_conservative_additive_attempt_failure_branch_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n154 = load_json(
        "fundamental_action_reconstruction/generated/n154_next_constructive_move_reduced_to_one_explicit_additive_success_failure_branch_split_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n154_binary_branch_split_explicit",
            "actual": n154["theorem_result"]["binary_branch_split_explicit"],
            "expected": True,
            "meaning": "N154 already freezes the binary verdict split",
        },
        {
            "id": "n154_failure_branch_name",
            "actual": n154["remaining_open_branches"][0],
            "expected": "future_failure_branch_discharge_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "the failure branch name is fixed",
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
        "stage": "F55",
        "lane": "conservative_failure_branch_first_only",
        "goal": "freeze_the_failure_branch_as_the_first_branch_to_attack_under_no_false_pass_ordering",
        "status": "F55_EXECUTED_FIRST_CONSERVATIVE_ADDITIVE_ATTEMPT_FAILURE_BRANCH_PACKET_NO_FALSE_PASS",
        "first_branch_to_attack": "explicit_failure_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
        "branch_ordering_basis": [
            "binary_split_already_frozen_by_N154",
            "no_constructed_source_object",
            "no_admissible_S_sel_int",
            "failure_first_is_methodologically_more_conservative_than_success_first",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F55",
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
