#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n128_last_positive_branch_reduced_to_one_initial_source_seed_target_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p117 = load_json(
        "fundamental_action_reconstruction/generated/p117_initial_future_strict_core_selector_source_seed_target_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p117_seed_target_reduced",
            "actual": p117["target_state"][
                "last_positive_branch_reduced_to_one_initial_seed_target"
            ],
            "expected": True,
            "meaning": "P117 already proves the last branch is reduced to one initial seed target",
        },
        {
            "id": "seed_source_name",
            "actual": p117["target_state"]["initial_seed_target"][
                "strict_core_source_object"
            ],
            "expected": "S_sel_int",
            "meaning": "the initial seed begins with the new source object",
        },
        {
            "id": "seed_orientation_name",
            "actual": p117["target_state"]["initial_seed_target"][
                "internal_orientation_export"
            ],
            "expected": "E_orient",
            "meaning": "the initial seed includes the orientation export",
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
            "step": "N128",
            "status": "N128_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_INITIAL_SOURCE_SEED_TARGET_STATE",
            "scope": "initial_future_source_seed_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N128",
            "status": "N128_DISCHARGED_LAST_POSITIVE_BRANCH_REDUCED_TO_ONE_INITIAL_SOURCE_SEED_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "initial_future_source_seed_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "initial_seed_target_reduced": True,
                "initial_seed_target": p117["target_state"]["initial_seed_target"],
                "downstream_chain_left_open": p117["target_state"][
                    "downstream_chain_left_open"
                ],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_construction_of_S_sel_int_and_E_orient_seed",
                "future_completion_of_B_sel_R_sel_O_sel_after_seed",
            ],
            "hard_limits": [
                "seed_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT)


if __name__ == "__main__":
    main()
