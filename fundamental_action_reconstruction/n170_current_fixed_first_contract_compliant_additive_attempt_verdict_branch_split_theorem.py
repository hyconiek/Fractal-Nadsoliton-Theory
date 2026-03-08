#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n170_current_fixed_first_contract_compliant_additive_attempt_verdict_branch_split_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p154 = load_json(
        "fundamental_action_reconstruction/generated/p154_first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p154_status",
            "actual": p154["status"],
            "expected": "CURRENT_REPO_REDUCES_THE_FIXED_FIRST_CONTRACT_COMPLIANT_FUTURE_ADDITIVE_UPSTREAM_SOURCE_CONSTRUCTION_ATTEMPT_VERDICT_TARGET_TO_ONE_EXPLICIT_SUCCESS_FAILURE_BRANCH_SPLIT_AFTER_P154",
            "meaning": "P154 already reduces the fixed verdict target to one explicit success/failure split",
        }
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
            "step": "N170",
            "status": "N170_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_VERDICT_BRANCH_SPLIT_STATE",
            "scope": "current_fixed_first_contract_compliant_additive_attempt_verdict_branch_split_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N170",
            "status": "N170_DISCHARGED_CURRENT_FIXED_FIRST_CONTRACT_COMPLIANT_ADDITIVE_ATTEMPT_VERDICT_BRANCH_SPLIT_THEOREM_NO_FALSE_PASS",
            "scope": "current_fixed_first_contract_compliant_additive_attempt_verdict_branch_split_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "verdict_branch_split_active": True,
                "allowed_branches": [
                    "explicit_success_verdict_for_construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
                    "explicit_failure_verdict_for_construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
                ],
                "full_closure_pass": False,
            },
            "hard_limits": [
                "no_success_verdict",
                "no_failure_verdict",
                "future_additive_source_object_not_yet_constructed",
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
