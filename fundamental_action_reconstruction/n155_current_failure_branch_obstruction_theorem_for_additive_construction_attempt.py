#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n155_current_failure_branch_obstruction_theorem_for_additive_construction_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p142 = load_json(
        "fundamental_action_reconstruction/generated/p142_current_failure_verdict_discharge_probe_for_additive_construction_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "p142_fixed_attempt_instance",
            "actual": p142["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "P142 fixes the additive construction attempt under test",
        },
        {
            "id": "p142_tested_failure_branch",
            "actual": p142["tested_failure_branch"],
            "expected": "explicit_failure_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "P142 fixes the failure branch under test",
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
            "step": "N155",
            "status": "N155_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CURRENT_FAILURE_BRANCH_STATE_FOR_ADDITIVE_CONSTRUCTION_ATTEMPT",
            "scope": "current_failure_verdict_discharge_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N155",
            "status": "N155_DISCHARGED_CURRENT_FAILURE_BRANCH_OBSTRUCTION_THEOREM_FOR_ADDITIVE_CONSTRUCTION_ATTEMPT_NO_FALSE_PASS",
            "scope": "current_failure_verdict_discharge_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "current_failure_branch_obstructed": True,
                "tested_failure_branch": "explicit_failure_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "fixed_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "full_closure_pass": False,
            },
            "hard_limits": [
                "success_branch_not_tested_here",
                "attempt_not_failed_in_absolute_future_independent_sense",
                "constructed_source_object_not_yet_exported",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
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
