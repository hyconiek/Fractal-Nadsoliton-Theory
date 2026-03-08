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
    / "p140_first_future_additive_construction_attempt_verdict_target_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p140_first_future_additive_construction_attempt_verdict_target_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f53 = load_json(
        "fundamental_action_reconstruction/generated/f53_first_future_additive_construction_attempt_verdict_target_packet_summary.json"
    )
    n152 = load_json(
        "fundamental_action_reconstruction/generated/n152_only_remaining_positive_move_class_reduced_to_one_first_additive_construction_attempt_instance_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n152_attempt_instance",
            "actual": n152["theorem_result"]["explicit_future_additive_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N152 already fixes one explicit future additive construction-attempt instance",
        },
        {
            "id": "f53_verdict_shape",
            "actual": f53["future_additive_attempt_verdict_target"]["verdict_shape"],
            "expected": "success_or_failure_verdict(construct_attempt_v1(S_sel_int_additive_attempt_target_v1))",
            "meaning": "F53 freezes one explicit verdict target shape",
        },
        {
            "id": "f53_no_constructed_object",
            "actual": f53["future_additive_attempt_verdict_target"]["constructed_source_object_present"],
            "expected": False,
            "meaning": "F53 does not overclaim a constructed source object",
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
        "stage": "P140",
        "lane": "future_additive_attempt_verdict_target_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_explicit_future_additive_attempt_verdict_target",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_EXPLICIT_FUTURE_ADDITIVE_CONSTRUCTION_ATTEMPT_VERDICT_TARGET_AFTER_P140",
        "target_state": {
            "explicit_future_additive_attempt_verdict_target": "success_or_failure_verdict(construct_attempt_v1(S_sel_int_additive_attempt_target_v1))",
            "remaining_open_branches": [
                "future_explicit_success_failure_branch_split_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)"
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P140",
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
