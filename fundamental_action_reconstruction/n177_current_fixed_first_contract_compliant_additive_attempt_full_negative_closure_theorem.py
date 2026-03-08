#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n177_current_fixed_first_contract_compliant_additive_attempt_full_negative_closure_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n171 = load_json(
        "fundamental_action_reconstruction/generated/n171_current_failure_branch_obstruction_theorem_for_fixed_first_contract_compliant_additive_attempt_summary.json"
    )
    n172 = load_json(
        "fundamental_action_reconstruction/generated/n172_current_success_branch_obstruction_theorem_for_fixed_first_contract_compliant_additive_attempt_summary.json"
    )
    n176 = load_json(
        "fundamental_action_reconstruction/generated/n176_current_fixed_first_contract_compliant_additive_attempt_post_verdict_lower_branch_full_negative_closure_theorem_summary.json"
    )

    expected_attempt = "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)"

    checks_spec = [
        {
            "id": "n171_failure_branch_obstructed",
            "actual": n171["theorem_result"]["failure_branch_current_state_obstructed"],
            "expected": True,
            "meaning": "N171 already obstructs the failure branch",
        },
        {
            "id": "n171_remaining_live_branch_name",
            "actual": n171["theorem_result"]["remaining_live_branch"],
            "expected": "explicit_success_verdict_for_construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "meaning": "N171 already fixes the remaining branch as success",
        },
        {
            "id": "n172_success_branch_obstructed",
            "actual": n172["theorem_result"]["success_branch_current_state_obstructed"],
            "expected": True,
            "meaning": "N172 already obstructs the success branch",
        },
        {
            "id": "n176_post_verdict_lower_closed",
            "actual": n176["theorem_result"]["contract_compliant_post_verdict_lower_branch_frontier_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N176 already closes the contract-compliant lower-branch frontier negatively",
        },
        {
            "id": "n176_fixed_attempt_instance",
            "actual": n176["theorem_result"]["fixed_attempt_instance"],
            "expected": expected_attempt,
            "meaning": "N176 already fixes the same contract-compliant additive attempt",
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
            "step": "N177",
            "status": "N177_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CONTRACT_COMPLIANT_FULL_NEGATIVE_CLOSURE_STATE",
            "scope": "current_fixed_first_contract_compliant_additive_attempt_full_negative_closure_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N177",
            "status": "N177_DISCHARGED_CURRENT_FIXED_FIRST_CONTRACT_COMPLIANT_ADDITIVE_ATTEMPT_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_fixed_first_contract_compliant_additive_attempt_full_negative_closure_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "fixed_first_contract_compliant_additive_attempt_closed_negatively_on_current_repo_state": True,
                "fixed_attempt_instance": expected_attempt,
                "full_closure_pass": False,
            },
            "hard_limits": [
                "constructed_source_object_not_yet_exported",
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
