#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n130_last_positive_branch_reduced_to_one_first_source_seed_construction_target_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p119 = load_json(
        "fundamental_action_reconstruction/generated/p119_first_source_seed_construction_target_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p119_first_source_seed_target_reduced",
            "actual": p119["target_state"][
                "last_positive_branch_reduced_to_one_first_source_seed_target"
            ],
            "expected": True,
            "meaning": "P119 already proves the last branch is reduced to one first source-seed construction target",
        },
        {
            "id": "source_seed_target_name",
            "actual": p119["target_state"]["first_source_seed_construction_target"][
                "strict_core_source_object"
            ],
            "expected": "S_sel_int",
            "meaning": "the first source-seed construction target is anchored on S_sel_int",
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
            "step": "N130",
            "status": "N130_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_SOURCE_SEED_TARGET_STATE",
            "scope": "first_source_seed_construction_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N130",
            "status": "N130_DISCHARGED_LAST_POSITIVE_BRANCH_REDUCED_TO_ONE_FIRST_SOURCE_SEED_CONSTRUCTION_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "first_source_seed_construction_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "last_positive_branch_reduced_to_one_first_source_seed_target": True,
                "first_source_seed_construction_target": p119["target_state"][
                    "first_source_seed_construction_target"
                ],
                "later_open_branches": p119["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_construction_of_admissible_S_sel_int",
                "future_derivation_of_admissible_E_orient_from_S_sel_int",
                "future_completion_of_B_sel_R_sel_O_sel_after_seed_package",
            ],
            "hard_limits": [
                "source_seed_not_yet_constructed",
                "orientation_export_not_yet_constructed",
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
