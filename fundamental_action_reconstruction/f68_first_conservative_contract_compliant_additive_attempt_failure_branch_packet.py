#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f68_first_conservative_contract_compliant_additive_attempt_failure_branch_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f67 = load_json(
        "fundamental_action_reconstruction/generated/f67_first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_refinement_packet_summary.json"
    )
    n170 = load_json(
        "fundamental_action_reconstruction/generated/n170_current_fixed_first_contract_compliant_additive_attempt_verdict_branch_split_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f67_branch_split_active",
            "actual": f67["branch_split"][
                "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_split_active"
            ],
            "expected": True,
            "meaning": "F67 already activates the explicit branch split",
        },
        {
            "id": "n170_branch_split_active",
            "actual": n170["theorem_result"]["verdict_branch_split_active"],
            "expected": True,
            "meaning": "N170 already fixes the branch split theorem-level",
        },
        {
            "id": "n170_has_failure_branch",
            "actual": "explicit_failure_verdict_for_construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)" in n170["theorem_result"]["allowed_branches"],
            "expected": True,
            "meaning": "N170 already includes the explicit failure branch",
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
            "stage": "F68",
            "lane": "first_conservative_contract_compliant_additive_attempt_failure_branch_only",
            "status": "F68_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FAILURE_BRANCH_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F68",
            "lane": "first_conservative_contract_compliant_additive_attempt_failure_branch_only",
            "goal": "freeze_the_conservative_failure_branch_for_the_fixed_first_contract_compliant_future_attempt",
            "status": "F68_EXECUTED_FIRST_CONSERVATIVE_CONTRACT_COMPLIANT_ADDITIVE_ATTEMPT_FAILURE_BRANCH_PACKET_NO_FALSE_PASS",
            "checks": checks,
            "selected_branch": {
                "first_conservative_contract_compliant_additive_attempt_failure_branch_active": True,
                "selected_branch": "explicit_failure_verdict_for_construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            },
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
