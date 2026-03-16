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
    / "f639_first_future_constructed_source_object_realization_verdict_target_packet_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f639_first_future_constructed_source_object_realization_verdict_target_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n529 = load_json(
        "fundamental_action_reconstruction/generated/n529_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_attempt_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    checks_spec = [
        {
            "id": "n529_realization_attempt_reduced",
            "actual": n529["theorem_result"][
                "next_constructive_move_reduced_to_one_first_future_realization_attempt"
            ],
            "expected": True,
            "meaning": "N529 reduces the next constructive move to one first future realization attempt (v1)",
        },
        {
            "id": "n529_attempt_name",
            "actual": n529["theorem_result"]["first_future_constructed_realization_attempt"][
                "attempt_name"
            ],
            "expected": "S_sel_int_new_object_constructed_realization_attempt_v1",
            "meaning": "the realization attempt instance name is fixed (v1)",
        },
        {
            "id": "n529_attempt_shape",
            "actual": n529["theorem_result"]["first_future_constructed_realization_attempt"][
                "attempt_shape"
            ],
            "expected": "realize_as_constructed_source_object_attempt_v1",
            "meaning": "the realization attempt shape is fixed (v1)",
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
        "target_name": "S_sel_int_new_object_constructed_realization_verdict_target_v1",
        "verdict_shape": "success_or_failure_verdict",
        "input_realization_attempt": "S_sel_int_new_object_constructed_realization_attempt_v1",
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
        "stage": "F639",
        "lane": "first_future_constructed_source_object_realization_verdict_target_only",
        "goal": "freeze_the_narrowest_first_verdict_target_above_the_fixed_seed_v1_realization_attempt_instance",
        "status": "F639_EXECUTED_FIRST_FUTURE_CONSTRUCTED_SOURCE_OBJECT_REALIZATION_VERDICT_TARGET_PACKET_FOR_SEED_V1_NO_FALSE_PASS",
        "first_future_realization_verdict_target": first_future_realization_verdict_target,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F639",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "first_future_realization_verdict_target": artifact[
            "first_future_realization_verdict_target"
        ],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

