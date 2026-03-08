#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f41_first_future_constructed_source_object_realization_target_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f41_first_future_constructed_source_object_realization_target_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n138 = load_json(
        "fundamental_action_reconstruction/generated/n138_next_constructive_move_reduced_to_one_first_future_genuinely_new_source_object_lift_bind_attempt_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n138_first_future_attempt_reduced",
            "actual": n138["theorem_result"][
                "next_constructive_move_reduced_to_one_first_future_attempt"
            ],
            "expected": True,
            "meaning": "N138 already reduces the next constructive move to one first future attempted construction instance",
        },
        {
            "id": "n138_attempt_name",
            "actual": n138["theorem_result"]["first_future_lift_bind_attempt"][
                "attempt_name"
            ],
            "expected": "S_sel_int_new_object_lift_bind_attempt_v0",
            "meaning": "the attempt instance name is fixed",
        },
        {
            "id": "n138_attempt_shape",
            "actual": n138["theorem_result"]["first_future_lift_bind_attempt"][
                "attempt_shape"
            ],
            "expected": "strict_core_single_object_lift_bind_attempt_v0",
            "meaning": "the attempt shape is already frozen",
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

    first_future_constructed_realization_target = {
        "target_name": "S_sel_int_new_object_constructed_realization_target_v0",
        "realization_shape": "realize_as_constructed_source_object",
        "input_attempt_instance": "S_sel_int_new_object_lift_bind_attempt_v0",
        "target_output": "future_constructed_source_object_for_S_sel_int",
        "counts_only_as": "first_future_constructed_source_object_realization_target",
        "does_not_count_as": [
            "constructed_source_object",
            "admissible_S_sel_int",
            "admissible_E_orient",
            "strict_core_selector_closure",
        ],
    }

    artifact = {
        "stage": "F41",
        "lane": "first_future_constructed_source_object_realization_target_only",
        "goal": "freeze_the_narrowest_first_realization_target_above_the_fixed_future_attempt_instance",
        "status": "F41_EXECUTED_FIRST_FUTURE_CONSTRUCTED_SOURCE_OBJECT_REALIZATION_TARGET_PACKET_NO_FALSE_PASS",
        "first_future_constructed_realization_target": first_future_constructed_realization_target,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F41",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "first_future_constructed_realization_target": artifact[
            "first_future_constructed_realization_target"
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
