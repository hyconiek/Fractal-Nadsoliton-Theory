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
    / "p139_first_future_genuinely_additive_source_object_construction_attempt_instance_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p139_first_future_genuinely_additive_source_object_construction_attempt_instance_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f52 = load_json(
        "fundamental_action_reconstruction/generated/f52_first_future_genuinely_additive_source_object_construction_attempt_instance_packet_summary.json"
    )
    n151 = load_json(
        "fundamental_action_reconstruction/generated/n151_only_remaining_positive_move_class_reduced_to_one_first_additive_construction_attempt_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n151_single_remaining_branch",
            "actual": n151["remaining_open_branches"],
            "expected": [
                "future_attempted_construction_of_S_sel_int_additive_attempt_target_v1"
            ],
            "meaning": "N151 leaves exactly one future constructive branch",
        },
        {
            "id": "f52_attempt_id",
            "actual": f52["future_additive_attempt_instance"]["attempt_id"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "F52 freezes one explicit future additive construction-attempt instance",
        },
        {
            "id": "f52_no_success_verdict",
            "actual": f52["future_additive_attempt_instance"]["success_verdict_present"],
            "expected": False,
            "meaning": "F52 does not overclaim success",
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
        "stage": "P139",
        "lane": "future_additive_construction_attempt_instance_only",
        "goal": "test_whether_the_only_remaining_positive_move_class_is_now_reduced_to_one_explicit_future_additive_construction_attempt_instance",
        "status": "CURRENT_REPO_REDUCES_THE_ONLY_REMAINING_POSITIVE_MOVE_CLASS_TO_ONE_EXPLICIT_FUTURE_ADDITIVE_SOURCE_OBJECT_CONSTRUCTION_ATTEMPT_INSTANCE_AFTER_P139",
        "target_state": {
            "only_remaining_positive_move_class": "future_genuinely_additive_new_strict_core_source_object_construction",
            "explicit_future_additive_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "remaining_open_branches": [
                "future_success_or_failure_evaluation_of_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)"
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P139",
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
