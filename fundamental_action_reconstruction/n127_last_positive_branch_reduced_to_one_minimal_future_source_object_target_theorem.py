#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n127_last_positive_branch_reduced_to_one_minimal_future_source_object_target_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p116 = load_json(
        "fundamental_action_reconstruction/generated/p116_future_strict_core_internal_selector_source_object_target_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p116_branch_reduced_to_one_target",
            "actual": p116["target_state"][
                "last_positive_branch_reduced_to_one_target"
            ],
            "expected": True,
            "meaning": "P116 already proves that the last positive branch is reduced to one target",
        },
        {
            "id": "target_source_object_name",
            "actual": p116["target_state"]["minimal_target_chain"][
                "strict_core_source_object"
            ],
            "expected": "S_sel_int",
            "meaning": "the reduced target starts with the new strict-core source object",
        },
        {
            "id": "target_operator_reachability_name",
            "actual": p116["target_state"]["minimal_target_chain"][
                "downstream_operator_reachability"
            ],
            "expected": "O_sel",
            "meaning": "the reduced target ends with downstream operator reachability",
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
            "step": "N127",
            "status": "N127_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FUTURE_SOURCE_OBJECT_TARGET_STATE",
            "scope": "last_positive_branch_reduction_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N127",
            "status": "N127_DISCHARGED_LAST_POSITIVE_BRANCH_REDUCED_TO_ONE_MINIMAL_FUTURE_SOURCE_OBJECT_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "last_positive_branch_reduction_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "last_positive_branch_reduced_to_one_target": True,
                "future_source_object_target_chain": p116["target_state"][
                    "minimal_target_chain"
                ],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_construction_of_S_sel_int_E_orient_B_sel_R_sel_O_sel_chain",
            ],
            "hard_limits": [
                "target_not_yet_constructed",
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
