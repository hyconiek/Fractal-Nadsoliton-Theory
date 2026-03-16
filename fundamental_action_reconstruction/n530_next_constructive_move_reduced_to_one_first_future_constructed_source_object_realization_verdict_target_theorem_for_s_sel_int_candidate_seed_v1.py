#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n530_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_verdict_target_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p639 = load_json(
        "fundamental_action_reconstruction/generated/p639_first_future_constructed_source_object_realization_verdict_target_probe_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    checks_spec = [
        {
            "id": "p639_first_future_verdict_target_reduced",
            "actual": p639["target_state"][
                "next_constructive_move_reduced_to_one_first_future_verdict_target"
            ],
            "expected": True,
            "meaning": "P639 proves the next constructive move is reduced to one first future verdict target (v1)",
        },
        {
            "id": "target_name",
            "actual": p639["target_state"]["first_future_realization_verdict_target"][
                "target_name"
            ],
            "expected": "S_sel_int_new_object_constructed_realization_verdict_target_v1",
            "meaning": "the first future realization-verdict target is explicitly named (v1)",
        },
        {
            "id": "verdict_shape",
            "actual": p639["target_state"]["first_future_realization_verdict_target"][
                "verdict_shape"
            ],
            "expected": "success_or_failure_verdict",
            "meaning": "the theorem stays scoped to the verdict-target stage rather than a verdict",
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
            "step": "N530",
            "status": "N530_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_FUTURE_VERDICT_TARGET_STATE_FOR_SEED_V1",
            "scope": "first_future_constructed_source_object_realization_verdict_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
    else:
        summary = {
            "step": "N530",
            "status": "N530_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_FIRST_FUTURE_CONSTRUCTED_SOURCE_OBJECT_REALIZATION_VERDICT_TARGET_THEOREM_FOR_SEED_V1_NO_FALSE_PASS",
            "scope": "first_future_constructed_source_object_realization_verdict_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_constructive_move_reduced_to_one_first_future_verdict_target": True,
                "first_future_realization_verdict_target": p639["target_state"][
                    "first_future_realization_verdict_target"
                ],
                "later_open_branches": p639["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": p639["target_state"]["later_open_branches"],
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

