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
    / "p635_future_genuinely_new_source_object_lift_bind_target_probe_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p635_future_genuinely_new_source_object_lift_bind_target_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n525 = load_json(
        "fundamental_action_reconstruction/generated/n525_current_first_clause_obstruction_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
    )
    f635 = load_json(
        "fundamental_action_reconstruction/generated/f635_future_genuinely_new_source_object_lift_bind_target_packet_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    reduced_to_one_future_lift_bind_target = (
        n525["theorem_result"]["first_clause_currently_satisfied"] is False
        and f635["future_lift_bind_target"]["future_target_name"]
        == "S_sel_int_new_object_target_v1"
        and f635["future_lift_bind_target"]["construction_shape"]
        == "strict_core_single_object_lift_bind"
    )

    checks_spec = [
        {
            "id": "n525_seed_v1_first_clause_blocked",
            "actual": n525["theorem_result"]["first_clause_currently_satisfied"],
            "expected": False,
            "meaning": "N525 blocks the seed-v1 attempt on the first clause",
        },
        {
            "id": "f635_future_target_name",
            "actual": f635["future_lift_bind_target"]["future_target_name"],
            "expected": "S_sel_int_new_object_target_v1",
            "meaning": "F635 exports one explicit future lift/bind target for seed v1",
        },
        {
            "id": "reduced_to_one_future_lift_bind_target",
            "actual": reduced_to_one_future_lift_bind_target,
            "expected": True,
            "meaning": "the next constructive move is reduced to one future lift/bind target (v1)",
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
        "stage": "P635",
        "lane": "future_genuinely_new_source_object_lift_bind_target_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_future_genuinely_new_source_object_lift_bind_target_for_seed_v1",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_FUTURE_GENUINELY_NEW_SOURCE_OBJECT_LIFT_BIND_TARGET_FOR_S_SEL_INT_CANDIDATE_SEED_V1_AFTER_P635",
        "target_state": {
            "next_constructive_move_reduced_to_one_future_lift_bind_target": reduced_to_one_future_lift_bind_target,
            "future_lift_bind_target": f635["future_lift_bind_target"],
            "later_open_branches": [
                "future_attempted_construction_of_S_sel_int_new_object_target_v1",
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
        "stage": "P635",
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

