#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f43_first_future_constructed_source_object_realization_verdict_target_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f43_first_future_constructed_source_object_realization_verdict_target_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n140 = load_json(
        "fundamental_action_reconstruction/generated/n140_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_attempt_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n140_realization_attempt_reduced",
            "actual": n140["theorem_result"][
                "next_constructive_move_reduced_to_one_first_future_realization_attempt"
            ],
            "expected": True,
            "meaning": "N140 already reduces the next constructive move to one first future realization attempt",
        },
        {
            "id": "n140_attempt_name",
            "actual": n140["theorem_result"]["first_future_constructed_realization_attempt"][
                "attempt_name"
            ],
            "expected": "S_sel_int_new_object_constructed_realization_attempt_v0",
            "meaning": "the realization attempt instance name is fixed",
        },
        {
            "id": "n140_attempt_shape",
            "actual": n140["theorem_result"]["first_future_constructed_realization_attempt"][
                "attempt_shape"
            ],
            "expected": "realize_as_constructed_source_object_attempt_v0",
            "meaning": "the realization attempt shape is already frozen",
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

    first_future_realization_verdict_target = {
        "target_name": "S_sel_int_new_object_constructed_realization_verdict_target_v0",
        "verdict_shape": "success_or_failure_verdict",
        "input_realization_attempt": "S_sel_int_new_object_constructed_realization_attempt_v0",
        "counts_only_as": "first_future_constructed_source_object_realization_verdict_target",
        "does_not_count_as": [
            "success_verdict",
            "failure_verdict",
            "constructed_source_object",
            "admissible_S_sel_int",
            "admissible_E_orient",
            "strict_core_selector_closure",
        ],
    }

    artifact = {
        "stage": "F43",
        "lane": "first_future_constructed_source_object_realization_verdict_target_only",
        "goal": "freeze_the_narrowest_first_verdict_target_above_the_fixed_future_realization_attempt_instance",
        "status": "F43_EXECUTED_FIRST_FUTURE_CONSTRUCTED_SOURCE_OBJECT_REALIZATION_VERDICT_TARGET_PACKET_NO_FALSE_PASS",
        "first_future_realization_verdict_target": first_future_realization_verdict_target,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F43",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "first_future_realization_verdict_target": artifact[
            "first_future_realization_verdict_target"
        ],
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
