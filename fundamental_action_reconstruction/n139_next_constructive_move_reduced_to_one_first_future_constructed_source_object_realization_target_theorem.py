#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n139_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_target_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p128 = load_json(
        "fundamental_action_reconstruction/generated/p128_first_future_constructed_source_object_realization_target_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p128_first_future_realization_target_reduced",
            "actual": p128["target_state"][
                "next_constructive_move_reduced_to_one_first_future_realization_target"
            ],
            "expected": True,
            "meaning": "P128 already proves the next constructive move is reduced to one first future realization target",
        },
        {
            "id": "target_name",
            "actual": p128["target_state"]["first_future_constructed_realization_target"][
                "target_name"
            ],
            "expected": "S_sel_int_new_object_constructed_realization_target_v0",
            "meaning": "the first future realization target is explicitly named",
        },
        {
            "id": "realization_shape",
            "actual": p128["target_state"]["first_future_constructed_realization_target"][
                "realization_shape"
            ],
            "expected": "realize_as_constructed_source_object",
            "meaning": "the theorem stays scoped to the realization-target stage rather than a realized object",
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
            "step": "N139",
            "status": "N139_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_FUTURE_REALIZATION_TARGET_STATE",
            "scope": "first_future_constructed_source_object_realization_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N139",
            "status": "N139_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_FIRST_FUTURE_CONSTRUCTED_SOURCE_OBJECT_REALIZATION_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "first_future_constructed_source_object_realization_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_constructive_move_reduced_to_one_first_future_realization_target": True,
                "first_future_constructed_realization_target": p128["target_state"][
                    "first_future_constructed_realization_target"
                ],
                "later_open_branches": p128["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_attempted_realization_of_S_sel_int_new_object_constructed_realization_target_v0_into_a_constructed_source_object",
                "future_admissibility_test_of_a_future_constructed_source_object",
                "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
                "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            ],
            "hard_limits": [
                "constructed_source_object_not_yet_realized",
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
