#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n178_current_fixed_first_contract_compliant_additive_attempt_nonreopening_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n149 = load_json(
        "fundamental_action_reconstruction/generated/n149_current_repo_constructive_selector_frontier_exhaustion_theorem_summary.json"
    )
    n166 = load_json(
        "fundamental_action_reconstruction/generated/n166_current_only_honest_positive_work_contract_theorem_summary.json"
    )
    n177 = load_json(
        "fundamental_action_reconstruction/generated/n177_current_fixed_first_contract_compliant_additive_attempt_full_negative_closure_theorem_summary.json"
    )

    expected_attempt = "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)"

    checks_spec = [
        {
            "id": "n149_constructive_selector_frontier_exhausted",
            "actual": n149["theorem_result"]["constructive_selector_frontier_exhausted_on_current_repo_state"],
            "expected": True,
            "meaning": "N149 already exhausts the constructive selector frontier on current exports",
        },
        {
            "id": "n166_only_honest_positive_work_contract_active",
            "actual": n166["theorem_result"]["only_honest_positive_work_contract_active"],
            "expected": True,
            "meaning": "N166 already restricts honest positive work to one explicit contract",
        },
        {
            "id": "n166_remaining_positive_work_must_be_genuinely_additive",
            "actual": n166["theorem_result"]["remaining_positive_work_must_be_genuinely_additive"],
            "expected": True,
            "meaning": "N166 already requires that remaining positive work be genuinely additive",
        },
        {
            "id": "n177_fixed_attempt_closed_negatively",
            "actual": n177["theorem_result"]["fixed_first_contract_compliant_additive_attempt_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N177 already closes the fixed first contract-compliant additive attempt negatively",
        },
        {
            "id": "n177_fixed_attempt_instance",
            "actual": n177["theorem_result"]["fixed_attempt_instance"],
            "expected": expected_attempt,
            "meaning": "N177 already fixes the same contract-compliant additive attempt",
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
            "step": "N178",
            "status": "N178_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NONREOPENING_STATE",
            "scope": "current_fixed_first_contract_compliant_additive_attempt_nonreopening_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N178",
            "status": "N178_DISCHARGED_CURRENT_FIXED_FIRST_CONTRACT_COMPLIANT_ADDITIVE_ATTEMPT_NONREOPENING_THEOREM_NO_FALSE_PASS",
            "scope": "current_fixed_first_contract_compliant_additive_attempt_nonreopening_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "current_fixed_first_contract_compliant_additive_attempt_does_not_reopen_constructive_selector_frontier": True,
                "constructive_selector_frontier_remains_exhausted_on_current_repo_state": True,
                "only_remaining_positive_work_still_future_and_additive": True,
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
