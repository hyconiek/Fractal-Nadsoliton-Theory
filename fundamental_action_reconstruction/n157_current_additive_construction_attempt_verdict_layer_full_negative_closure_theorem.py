#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n157_current_additive_construction_attempt_verdict_layer_full_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n155 = load_json(
        "fundamental_action_reconstruction/generated/n155_current_failure_branch_obstruction_theorem_for_additive_construction_attempt_summary.json"
    )
    n156 = load_json(
        "fundamental_action_reconstruction/generated/n156_current_success_branch_obstruction_theorem_for_additive_construction_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "n155_failure_branch_closed",
            "actual": n155["theorem_result"]["current_failure_branch_obstructed"],
            "expected": True,
            "meaning": "N155 already closes the failure branch negatively on the current repo state",
        },
        {
            "id": "n156_success_branch_closed",
            "actual": n156["theorem_result"]["current_success_branch_obstructed"],
            "expected": True,
            "meaning": "N156 already closes the success branch negatively on the current repo state",
        },
        {
            "id": "n155_fixed_attempt_instance",
            "actual": n155["theorem_result"]["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N155 fixes the additive construction attempt under test",
        },
        {
            "id": "n156_fixed_attempt_instance",
            "actual": n156["theorem_result"]["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N156 fixes the same additive construction attempt under test",
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
            "step": "N157",
            "status": "N157_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADDITIVE_VERDICT_LAYER_STATE",
            "scope": "current_additive_construction_attempt_verdict_layer_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N157",
            "status": "N157_DISCHARGED_CURRENT_ADDITIVE_CONSTRUCTION_ATTEMPT_VERDICT_LAYER_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_additive_construction_attempt_verdict_layer_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "fixed_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "failure_branch_closed_negatively_on_current_repo_state": True,
                "success_branch_closed_negatively_on_current_repo_state": True,
                "additive_construction_attempt_verdict_layer_closed_negatively_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [],
            "hard_limits": [
                "no_proof_that_future_additive_source_object_construction_is_impossible_forever",
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
