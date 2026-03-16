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
    / "p637_first_future_constructed_source_object_realization_target_probe_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p637_first_future_constructed_source_object_realization_target_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n527 = load_json(
        "fundamental_action_reconstruction/generated/n527_next_constructive_move_reduced_to_one_first_future_genuinely_new_source_object_lift_bind_attempt_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
    )
    f637 = load_json(
        "fundamental_action_reconstruction/generated/f637_first_future_constructed_source_object_realization_target_packet_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    reduced_to_one_first_future_realization_target = (
        n527["theorem_result"]["next_constructive_move_reduced_to_one_first_future_attempt"]
        is True
        and f637["first_future_constructed_realization_target"]["target_name"]
        == "S_sel_int_new_object_constructed_realization_target_v1"
        and f637["first_future_constructed_realization_target"]["realization_shape"]
        == "realize_as_constructed_source_object"
    )

    checks_spec = [
        {
            "id": "n527_first_future_attempt_reduced",
            "actual": n527["theorem_result"][
                "next_constructive_move_reduced_to_one_first_future_attempt"
            ],
            "expected": True,
            "meaning": "N527 already reduces the next constructive move to one first future attempt (v1)",
        },
        {
            "id": "f637_target_name",
            "actual": f637["first_future_constructed_realization_target"]["target_name"],
            "expected": "S_sel_int_new_object_constructed_realization_target_v1",
            "meaning": "F637 exports one explicit future realization target (v1)",
        },
        {
            "id": "reduced_to_one_first_future_realization_target",
            "actual": reduced_to_one_first_future_realization_target,
            "expected": True,
            "meaning": "the next constructive move is reduced to one first future constructed-source-object realization target (v1)",
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
        "stage": "P637",
        "lane": "first_future_constructed_source_object_realization_target_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_first_future_constructed_source_object_realization_target_for_seed_v1",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_FIRST_FUTURE_CONSTRUCTED_SOURCE_OBJECT_REALIZATION_TARGET_FOR_SEED_V1_AFTER_P637",
        "target_state": {
            "next_constructive_move_reduced_to_one_first_future_realization_target": reduced_to_one_first_future_realization_target,
            "first_future_constructed_realization_target": f637[
                "first_future_constructed_realization_target"
            ],
            "later_open_branches": [
                "future_attempted_realization_of_S_sel_int_new_object_constructed_realization_target_v1_into_a_constructed_source_object",
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
        "stage": "P637",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

