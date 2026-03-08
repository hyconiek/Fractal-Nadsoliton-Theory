#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n172_current_success_branch_obstruction_theorem_for_fixed_first_contract_compliant_additive_attempt_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p156 = load_json(
        "fundamental_action_reconstruction/generated/p156_current_success_verdict_discharge_probe_for_fixed_first_contract_compliant_additive_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "p156_status",
            "actual": p156["status"],
            "expected": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_SUCCESS_VERDICT_DISCHARGE_FOR_THE_FIXED_FIRST_CONTRACT_COMPLIANT_FUTURE_ADDITIVE_UPSTREAM_SOURCE_CONSTRUCTION_ATTEMPT_AFTER_P156",
            "meaning": "P156 already keeps the success verdict branch obstructed",
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
            "step": "N172",
            "status": "N172_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SUCCESS_BRANCH_OBSTRUCTION_STATE",
            "scope": "current_success_branch_obstruction_for_fixed_first_contract_compliant_additive_attempt_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N172",
            "status": "N172_DISCHARGED_CURRENT_SUCCESS_BRANCH_OBSTRUCTION_THEOREM_FOR_FIXED_FIRST_CONTRACT_COMPLIANT_ADDITIVE_ATTEMPT_NO_FALSE_PASS",
            "scope": "current_success_branch_obstruction_for_fixed_first_contract_compliant_additive_attempt_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "success_branch_current_state_obstructed": True,
                "binary_verdict_layer_fully_negative_on_current_repo_state": True,
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
