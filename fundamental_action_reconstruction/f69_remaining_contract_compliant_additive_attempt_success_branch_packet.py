#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f69_remaining_contract_compliant_additive_attempt_success_branch_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n171 = load_json(
        "fundamental_action_reconstruction/generated/n171_current_failure_branch_obstruction_theorem_for_fixed_first_contract_compliant_additive_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "n171_failure_branch_obstructed",
            "actual": n171["theorem_result"]["failure_branch_current_state_obstructed"],
            "expected": True,
            "meaning": "N171 already obstructs the fixed failure branch",
        },
        {
            "id": "n171_remaining_live_branch",
            "actual": n171["theorem_result"]["remaining_live_branch"],
            "expected": "explicit_success_verdict_for_construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "meaning": "N171 already fixes the remaining live branch as success",
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
            "stage": "F69",
            "status": "F69_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SUCCESS_BRANCH_REMAINING_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F69",
            "status": "F69_EXECUTED_REMAINING_CONTRACT_COMPLIANT_ADDITIVE_ATTEMPT_SUCCESS_BRANCH_PACKET_NO_FALSE_PASS",
            "selected_branch": {
                "remaining_contract_compliant_additive_attempt_success_branch_active": True,
                "selected_branch": "explicit_success_verdict_for_construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            },
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
