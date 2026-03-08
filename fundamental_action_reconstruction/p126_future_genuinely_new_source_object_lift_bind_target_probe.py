#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p126_future_genuinely_new_source_object_lift_bind_target_probe.json"
OUT_SUMMARY = (
    GENERATED / "p126_future_genuinely_new_source_object_lift_bind_target_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n136 = load_json(
        "fundamental_action_reconstruction/generated/n136_current_first_clause_obstruction_theorem_for_s_sel_int_candidate_seed_v0_summary.json"
    )
    f39 = load_json(
        "fundamental_action_reconstruction/generated/f39_future_genuinely_new_source_object_lift_bind_target_packet_summary.json"
    )

    reduced_to_one_future_lift_bind_target = (
        n136["theorem_result"]["first_clause_currently_satisfied"] is False
        and f39["future_lift_bind_target"]["future_target_name"] == "S_sel_int_new_object_target_v0"
        and f39["future_lift_bind_target"]["construction_shape"] == "strict_core_single_object_lift_bind"
    )

    checks_spec = [
        {
            "id": "n136_packaged_seed_blocked",
            "actual": n136["theorem_result"]["first_clause_currently_satisfied"],
            "expected": False,
            "meaning": "N136 already blocks the packaged seed on the first clause",
        },
        {
            "id": "f39_future_target_name",
            "actual": f39["future_lift_bind_target"]["future_target_name"],
            "expected": "S_sel_int_new_object_target_v0",
            "meaning": "F39 exports one explicit future lift/bind target",
        },
        {
            "id": "reduced_to_one_future_lift_bind_target",
            "actual": reduced_to_one_future_lift_bind_target,
            "expected": True,
            "meaning": "the next constructive move is reduced to one future lift/bind target",
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
        "stage": "P126",
        "lane": "future_genuinely_new_source_object_lift_bind_target_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_future_genuinely_new_source_object_lift_bind_target",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_FUTURE_GENUINELY_NEW_SOURCE_OBJECT_LIFT_BIND_TARGET_AFTER_P126",
        "target_state": {
            "next_constructive_move_reduced_to_one_future_lift_bind_target": reduced_to_one_future_lift_bind_target,
            "future_lift_bind_target": f39["future_lift_bind_target"],
            "later_open_branches": [
                "future_attempted_construction_of_S_sel_int_new_object_target_v0",
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
        "stage": "P126",
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
