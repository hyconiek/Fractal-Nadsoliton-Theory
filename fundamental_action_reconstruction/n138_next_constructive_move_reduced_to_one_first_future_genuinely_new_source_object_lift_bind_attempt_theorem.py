#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n138_next_constructive_move_reduced_to_one_first_future_genuinely_new_source_object_lift_bind_attempt_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p127 = load_json(
        "fundamental_action_reconstruction/generated/p127_first_future_genuinely_new_source_object_lift_bind_attempt_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p127_first_future_attempt_reduced",
            "actual": p127["target_state"][
                "next_constructive_move_reduced_to_one_first_future_attempt"
            ],
            "expected": True,
            "meaning": "P127 already proves the next constructive move is reduced to one first future attempt",
        },
        {
            "id": "attempt_name",
            "actual": p127["target_state"]["first_future_lift_bind_attempt"][
                "attempt_name"
            ],
            "expected": "S_sel_int_new_object_lift_bind_attempt_v0",
            "meaning": "the first future attempted construction instance is explicitly named",
        },
        {
            "id": "attempt_shape",
            "actual": p127["target_state"]["first_future_lift_bind_attempt"][
                "attempt_shape"
            ],
            "expected": "strict_core_single_object_lift_bind_attempt_v0",
            "meaning": "the theorem stays scoped to the attempt instance rather than a realized object",
        },
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "step": "N138",
            "status": "N138_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_FUTURE_ATTEMPT_STATE",
            "scope": "first_future_genuinely_new_source_object_lift_bind_attempt_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N138",
            "status": "N138_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_FIRST_FUTURE_GENUINELY_NEW_SOURCE_OBJECT_LIFT_BIND_ATTEMPT_THEOREM_NO_FALSE_PASS",
            "scope": "first_future_genuinely_new_source_object_lift_bind_attempt_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_constructive_move_reduced_to_one_first_future_attempt": True,
                "first_future_lift_bind_attempt": p127["target_state"][
                    "first_future_lift_bind_attempt"
                ],
                "later_open_branches": p127["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_attempted_realization_of_S_sel_int_new_object_lift_bind_attempt_v0_as_constructed_source_object",
                "future_admissibility_test_of_a_future_constructed_new_source_object",
                "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
                "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            ],
            "hard_limits": [
                "attempt_instance_not_yet_realized_as_constructed_object",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
