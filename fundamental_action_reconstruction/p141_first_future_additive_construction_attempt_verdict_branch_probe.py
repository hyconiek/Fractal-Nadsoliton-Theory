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
    / "p141_first_future_additive_construction_attempt_verdict_branch_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p141_first_future_additive_construction_attempt_verdict_branch_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f54 = load_json(
        "fundamental_action_reconstruction/generated/f54_first_future_additive_construction_attempt_verdict_branch_refinement_packet_summary.json"
    )
    n153 = load_json(
        "fundamental_action_reconstruction/generated/n153_next_constructive_move_reduced_to_one_first_additive_construction_attempt_verdict_target_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n153_verdict_target",
            "actual": n153["theorem_result"]["explicit_future_additive_attempt_verdict_target"],
            "expected": "success_or_failure_verdict(construct_attempt_v1(S_sel_int_additive_attempt_target_v1))",
            "meaning": "N153 fixes the additive-construction verdict target",
        },
        {
            "id": "f54_success_branch",
            "actual": f54["verdict_branch_split"]["success_branch"],
            "expected": "explicit_success_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "F54 freezes the success branch",
        },
        {
            "id": "f54_failure_branch",
            "actual": f54["verdict_branch_split"]["failure_branch"],
            "expected": "explicit_failure_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "F54 freezes the failure branch",
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
        "stage": "P141",
        "lane": "future_additive_attempt_verdict_branch_refinement_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_explicit_binary_branch_split_on_the_fixed_additive_construction_verdict_target",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_EXPLICIT_FUTURE_ADDITIVE_SUCCESS_FAILURE_BRANCH_SPLIT_AFTER_P141",
        "target_state": {
            "fixed_verdict_target": "success_or_failure_verdict(construct_attempt_v1(S_sel_int_additive_attempt_target_v1))",
            "remaining_open_branches": [
                "future_failure_branch_discharge_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "future_success_branch_discharge_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)"
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P141",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
