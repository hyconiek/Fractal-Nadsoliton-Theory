#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n156_current_success_branch_obstruction_theorem_for_additive_construction_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p143 = load_json(
        "fundamental_action_reconstruction/generated/p143_current_success_verdict_discharge_probe_for_additive_construction_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "p143_fixed_attempt_instance",
            "actual": p143["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "P143 fixes the additive construction attempt under test",
        },
        {
            "id": "p143_tested_success_branch",
            "actual": p143["tested_success_branch"],
            "expected": "explicit_success_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "P143 fixes the success branch under test",
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
            "step": "N156",
            "status": "N156_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CURRENT_SUCCESS_BRANCH_STATE_FOR_ADDITIVE_CONSTRUCTION_ATTEMPT",
            "scope": "current_success_verdict_discharge_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N156",
            "status": "N156_DISCHARGED_CURRENT_SUCCESS_BRANCH_OBSTRUCTION_THEOREM_FOR_ADDITIVE_CONSTRUCTION_ATTEMPT_NO_FALSE_PASS",
            "scope": "current_success_verdict_discharge_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "current_success_branch_obstructed": True,
                "tested_success_branch": "explicit_success_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "fixed_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "full_closure_pass": False,
            },
            "hard_limits": [
                "failure_branch_not_future_independent_impossibility_here",
                "attempt_not_succeeded",
                "attempt_not_failed",
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
