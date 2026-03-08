#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f40_first_future_genuinely_new_source_object_lift_bind_attempt_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f40_first_future_genuinely_new_source_object_lift_bind_attempt_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n137 = load_json(
        "fundamental_action_reconstruction/generated/n137_next_constructive_move_reduced_to_one_future_genuinely_new_source_object_lift_bind_target_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n137_future_target_reduced",
            "actual": n137["theorem_result"][
                "next_constructive_move_reduced_to_one_future_lift_bind_target"
            ],
            "expected": True,
            "meaning": "N137 already reduces the next constructive move to one future lift/bind target",
        },
        {
            "id": "n137_target_name",
            "actual": n137["theorem_result"]["future_lift_bind_target"][
                "future_target_name"
            ],
            "expected": "S_sel_int_new_object_target_v0",
            "meaning": "the future target name is fixed",
        },
        {
            "id": "n137_construction_shape",
            "actual": n137["theorem_result"]["future_lift_bind_target"][
                "construction_shape"
            ],
            "expected": "strict_core_single_object_lift_bind",
            "meaning": "the target shape is already frozen as a single-object lift/bind target",
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

    first_future_lift_bind_attempt = {
        "attempt_name": "S_sel_int_new_object_lift_bind_attempt_v0",
        "attempt_shape": "strict_core_single_object_lift_bind_attempt_v0",
        "input_materials": [
            "QW-2206_local_topological_protection_layer",
            "sigma_int_candidate",
        ],
        "intended_target_output": "future_S_sel_int",
        "counts_only_as": "first_future_genuinely_new_source_object_lift_bind_attempt",
        "does_not_count_as": [
            "constructed_source_object",
            "admissible_S_sel_int",
            "admissible_E_orient",
            "strict_core_selector_closure",
        ],
    }

    artifact = {
        "stage": "F40",
        "lane": "first_future_genuinely_new_source_object_lift_bind_attempt_only",
        "goal": "freeze_the_narrowest_first_attempted_construction_instance_on_the_fixed_future_genuinely_new_source_object_lift_bind_target",
        "status": "F40_EXECUTED_FIRST_FUTURE_GENUINELY_NEW_SOURCE_OBJECT_LIFT_BIND_ATTEMPT_PACKET_NO_FALSE_PASS",
        "first_future_lift_bind_attempt": first_future_lift_bind_attempt,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F40",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "first_future_lift_bind_attempt": artifact["first_future_lift_bind_attempt"],
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
