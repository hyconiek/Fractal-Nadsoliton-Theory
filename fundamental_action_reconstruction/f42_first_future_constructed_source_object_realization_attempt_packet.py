#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f42_first_future_constructed_source_object_realization_attempt_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f42_first_future_constructed_source_object_realization_attempt_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n139 = load_json(
        "fundamental_action_reconstruction/generated/n139_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_target_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n139_realization_target_reduced",
            "actual": n139["theorem_result"][
                "next_constructive_move_reduced_to_one_first_future_realization_target"
            ],
            "expected": True,
            "meaning": "N139 already reduces the next constructive move to one first future realization target",
        },
        {
            "id": "n139_target_name",
            "actual": n139["theorem_result"]["first_future_constructed_realization_target"][
                "target_name"
            ],
            "expected": "S_sel_int_new_object_constructed_realization_target_v0",
            "meaning": "the realization target name is fixed",
        },
        {
            "id": "n139_input_attempt_instance",
            "actual": n139["theorem_result"]["first_future_constructed_realization_target"][
                "input_attempt_instance"
            ],
            "expected": "S_sel_int_new_object_lift_bind_attempt_v0",
            "meaning": "the input attempt instance is already frozen",
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

    first_future_constructed_realization_attempt = {
        "attempt_name": "S_sel_int_new_object_constructed_realization_attempt_v0",
        "attempt_shape": "realize_as_constructed_source_object_attempt_v0",
        "input_attempt_instance": "S_sel_int_new_object_lift_bind_attempt_v0",
        "intended_target_output": "future_constructed_source_object_for_S_sel_int",
        "counts_only_as": "first_future_constructed_source_object_realization_attempt",
        "does_not_count_as": [
            "constructed_source_object",
            "admissible_S_sel_int",
            "admissible_E_orient",
            "strict_core_selector_closure",
        ],
    }

    artifact = {
        "stage": "F42",
        "lane": "first_future_constructed_source_object_realization_attempt_only",
        "goal": "freeze_the_narrowest_first_attempted_realization_instance_on_the_fixed_future_constructed_source_object_realization_target",
        "status": "F42_EXECUTED_FIRST_FUTURE_CONSTRUCTED_SOURCE_OBJECT_REALIZATION_ATTEMPT_PACKET_NO_FALSE_PASS",
        "first_future_constructed_realization_attempt": first_future_constructed_realization_attempt,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F42",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "first_future_constructed_realization_attempt": artifact[
            "first_future_constructed_realization_attempt"
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
