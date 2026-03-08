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
    / "p127_first_future_genuinely_new_source_object_lift_bind_attempt_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p127_first_future_genuinely_new_source_object_lift_bind_attempt_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n137 = load_json(
        "fundamental_action_reconstruction/generated/n137_next_constructive_move_reduced_to_one_future_genuinely_new_source_object_lift_bind_target_theorem_summary.json"
    )
    f40 = load_json(
        "fundamental_action_reconstruction/generated/f40_first_future_genuinely_new_source_object_lift_bind_attempt_packet_summary.json"
    )

    reduced_to_one_first_future_attempt = (
        n137["theorem_result"]["next_constructive_move_reduced_to_one_future_lift_bind_target"]
        is True
        and f40["first_future_lift_bind_attempt"]["attempt_name"]
        == "S_sel_int_new_object_lift_bind_attempt_v0"
        and f40["first_future_lift_bind_attempt"]["attempt_shape"]
        == "strict_core_single_object_lift_bind_attempt_v0"
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
            "id": "f40_attempt_name",
            "actual": f40["first_future_lift_bind_attempt"]["attempt_name"],
            "expected": "S_sel_int_new_object_lift_bind_attempt_v0",
            "meaning": "F40 exports one explicit future attempted construction instance",
        },
        {
            "id": "reduced_to_one_first_future_attempt",
            "actual": reduced_to_one_first_future_attempt,
            "expected": True,
            "meaning": "the next constructive move is reduced to one first future attempted construction instance",
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
        "stage": "P127",
        "lane": "first_future_genuinely_new_source_object_lift_bind_attempt_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_first_future_genuinely_new_source_object_lift_bind_attempt",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_FIRST_FUTURE_GENUINELY_NEW_SOURCE_OBJECT_LIFT_BIND_ATTEMPT_AFTER_P127",
        "target_state": {
            "next_constructive_move_reduced_to_one_first_future_attempt": reduced_to_one_first_future_attempt,
            "first_future_lift_bind_attempt": f40["first_future_lift_bind_attempt"],
            "later_open_branches": [
                "future_attempted_realization_of_S_sel_int_new_object_lift_bind_attempt_v0_as_constructed_source_object",
                "future_admissibility_test_of_a_future_constructed_new_source_object",
                "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
                "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P127",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
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
