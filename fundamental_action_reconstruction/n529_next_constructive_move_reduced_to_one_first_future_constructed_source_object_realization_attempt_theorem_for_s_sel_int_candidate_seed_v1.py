#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n529_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_attempt_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p638 = load_json(
        "fundamental_action_reconstruction/generated/p638_first_future_constructed_source_object_realization_attempt_probe_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    checks_spec = [
        {
            "id": "p638_first_future_realization_attempt_reduced",
            "actual": p638["target_state"][
                "next_constructive_move_reduced_to_one_first_future_realization_attempt"
            ],
            "expected": True,
            "meaning": "P638 proves the next constructive move is reduced to one first future realization attempt (v1)",
        },
        {
            "id": "attempt_name",
            "actual": p638["target_state"]["first_future_constructed_realization_attempt"][
                "attempt_name"
            ],
            "expected": "S_sel_int_new_object_constructed_realization_attempt_v1",
            "meaning": "the first future realization attempt instance is explicitly named (v1)",
        },
        {
            "id": "attempt_shape",
            "actual": p638["target_state"]["first_future_constructed_realization_attempt"][
                "attempt_shape"
            ],
            "expected": "realize_as_constructed_source_object_attempt_v1",
            "meaning": "the theorem stays scoped to the attempt stage rather than a realized object",
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
            "step": "N529",
            "status": "N529_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_FUTURE_REALIZATION_ATTEMPT_STATE_FOR_SEED_V1",
            "scope": "first_future_constructed_source_object_realization_attempt_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
    else:
        summary = {
            "step": "N529",
            "status": "N529_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_FIRST_FUTURE_CONSTRUCTED_SOURCE_OBJECT_REALIZATION_ATTEMPT_THEOREM_FOR_SEED_V1_NO_FALSE_PASS",
            "scope": "first_future_constructed_source_object_realization_attempt_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_constructive_move_reduced_to_one_first_future_realization_attempt": True,
                "first_future_constructed_realization_attempt": p638["target_state"][
                    "first_future_constructed_realization_attempt"
                ],
                "later_open_branches": p638["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": p638["target_state"]["later_open_branches"],
            "hard_limits": [
                "success_or_failure_verdict_not_yet_discharged",
                "constructed_source_object_not_yet_exported",
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

