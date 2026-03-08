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
    / "p129_first_future_constructed_source_object_realization_attempt_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p129_first_future_constructed_source_object_realization_attempt_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n139 = load_json(
        "fundamental_action_reconstruction/generated/n139_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_target_theorem_summary.json"
    )
    f42 = load_json(
        "fundamental_action_reconstruction/generated/f42_first_future_constructed_source_object_realization_attempt_packet_summary.json"
    )

    reduced_to_one_first_future_realization_attempt = (
        n139["theorem_result"][
            "next_constructive_move_reduced_to_one_first_future_realization_target"
        ]
        is True
        and f42["first_future_constructed_realization_attempt"]["attempt_name"]
        == "S_sel_int_new_object_constructed_realization_attempt_v0"
        and f42["first_future_constructed_realization_attempt"]["attempt_shape"]
        == "realize_as_constructed_source_object_attempt_v0"
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
            "id": "f42_attempt_name",
            "actual": f42["first_future_constructed_realization_attempt"][
                "attempt_name"
            ],
            "expected": "S_sel_int_new_object_constructed_realization_attempt_v0",
            "meaning": "F42 exports one explicit future realization attempt instance",
        },
        {
            "id": "reduced_to_one_first_future_realization_attempt",
            "actual": reduced_to_one_first_future_realization_attempt,
            "expected": True,
            "meaning": "the next constructive move is reduced to one first future realization attempt instance",
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
        "stage": "P129",
        "lane": "first_future_constructed_source_object_realization_attempt_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_first_future_constructed_source_object_realization_attempt",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_FIRST_FUTURE_CONSTRUCTED_SOURCE_OBJECT_REALIZATION_ATTEMPT_AFTER_P129",
        "target_state": {
            "next_constructive_move_reduced_to_one_first_future_realization_attempt": reduced_to_one_first_future_realization_attempt,
            "first_future_constructed_realization_attempt": f42[
                "first_future_constructed_realization_attempt"
            ],
            "later_open_branches": [
                "future_success_or_failure_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0",
                "future_admissibility_test_of_a_future_constructed_source_object",
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
        "stage": "P129",
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
