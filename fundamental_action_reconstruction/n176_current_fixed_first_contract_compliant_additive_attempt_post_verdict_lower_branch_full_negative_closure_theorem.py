#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n176_current_fixed_first_contract_compliant_additive_attempt_post_verdict_lower_branch_full_negative_closure_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n173 = load_json(
        "fundamental_action_reconstruction/generated/n173_current_contract_compliant_post_verdict_admissibility_branch_obstruction_theorem_summary.json"
    )
    n174 = load_json(
        "fundamental_action_reconstruction/generated/n174_current_contract_compliant_orientation_export_branch_obstruction_theorem_summary.json"
    )
    n175 = load_json(
        "fundamental_action_reconstruction/generated/n175_current_contract_compliant_downstream_completion_branch_obstruction_theorem_summary.json"
    )

    expected_attempt = "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)"

    checks_spec = [
        {
            "id": "n173_fixed_attempt_instance",
            "actual": n173["theorem_result"]["fixed_attempt_instance"],
            "expected": expected_attempt,
            "meaning": "N173 already fixes the contract-compliant additive attempt",
        },
        {
            "id": "n174_fixed_attempt_instance",
            "actual": n174["theorem_result"]["fixed_attempt_instance"],
            "expected": expected_attempt,
            "meaning": "N174 already fixes the same contract-compliant additive attempt",
        },
        {
            "id": "n175_fixed_attempt_instance",
            "actual": n175["theorem_result"]["fixed_attempt_instance"],
            "expected": expected_attempt,
            "meaning": "N175 already fixes the same contract-compliant additive attempt",
        },
        {
            "id": "n173_admissibility_obstructed",
            "actual": n173["theorem_result"]["current_contract_compliant_post_verdict_admissibility_branch_obstructed"],
            "expected": True,
            "meaning": "N173 already obstructs the contract-compliant admissibility branch",
        },
        {
            "id": "n174_orientation_obstructed",
            "actual": n174["theorem_result"]["current_contract_compliant_orientation_export_branch_obstructed"],
            "expected": True,
            "meaning": "N174 already obstructs the contract-compliant orientation-export branch",
        },
        {
            "id": "n175_downstream_obstructed",
            "actual": n175["theorem_result"]["current_contract_compliant_downstream_completion_branch_obstructed"],
            "expected": True,
            "meaning": "N175 already obstructs the contract-compliant downstream-completion branch",
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
            "step": "N176",
            "status": "N176_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CONTRACT_COMPLIANT_POST_VERDICT_LOWER_BRANCH_STATE",
            "scope": "current_fixed_first_contract_compliant_additive_attempt_post_verdict_lower_branch_full_negative_closure_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N176",
            "status": "N176_DISCHARGED_CURRENT_FIXED_FIRST_CONTRACT_COMPLIANT_ADDITIVE_ATTEMPT_POST_VERDICT_LOWER_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_fixed_first_contract_compliant_additive_attempt_post_verdict_lower_branch_full_negative_closure_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "contract_compliant_post_verdict_lower_branch_frontier_closed_negatively_on_current_repo_state": True,
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
