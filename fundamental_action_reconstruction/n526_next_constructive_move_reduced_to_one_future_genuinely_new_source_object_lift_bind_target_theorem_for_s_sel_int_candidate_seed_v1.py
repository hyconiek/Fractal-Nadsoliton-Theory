#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n526_next_constructive_move_reduced_to_one_future_genuinely_new_source_object_lift_bind_target_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p635 = load_json(
        "fundamental_action_reconstruction/generated/p635_future_genuinely_new_source_object_lift_bind_target_probe_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    checks_spec = [
        {
            "id": "p635_future_target_reduced",
            "actual": p635["target_state"][
                "next_constructive_move_reduced_to_one_future_lift_bind_target"
            ],
            "expected": True,
            "meaning": "P635 reduces the next constructive move to one explicit future lift/bind target",
        },
        {
            "id": "future_target_name",
            "actual": p635["target_state"]["future_lift_bind_target"]["future_target_name"],
            "expected": "S_sel_int_new_object_target_v1",
            "meaning": "the future target is explicitly named",
        },
        {
            "id": "construction_shape",
            "actual": p635["target_state"]["future_lift_bind_target"]["construction_shape"],
            "expected": "strict_core_single_object_lift_bind",
            "meaning": "the theorem stays scoped to a single-object lift/bind target",
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
            "step": "N526",
            "status": "N526_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FUTURE_LIFT_BIND_TARGET_STATE",
            "scope": "future_genuinely_new_source_object_lift_bind_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
    else:
        summary = {
            "step": "N526",
            "status": "N526_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_FUTURE_GENUINELY_NEW_SOURCE_OBJECT_LIFT_BIND_TARGET_THEOREM_FOR_SEED_V1_NO_FALSE_PASS",
            "scope": "future_genuinely_new_source_object_lift_bind_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_constructive_move_reduced_to_one_future_lift_bind_target": True,
                "future_lift_bind_target": p635["target_state"]["future_lift_bind_target"],
                "later_open_branches": p635["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_attempted_construction_of_S_sel_int_new_object_target_v1",
                "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
                "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            ],
            "hard_limits": [
                "future_lift_bind_target_not_yet_constructed",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

