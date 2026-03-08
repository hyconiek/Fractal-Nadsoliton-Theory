#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p155_current_failure_verdict_discharge_probe_for_fixed_first_contract_compliant_additive_attempt_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f68 = load_json(
        "fundamental_action_reconstruction/generated/f68_first_conservative_contract_compliant_additive_attempt_failure_branch_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "f68_failure_branch_active",
            "actual": f68["selected_branch"][
                "first_conservative_contract_compliant_additive_attempt_failure_branch_active"
            ],
            "expected": True,
            "meaning": "F68 already activates the conservative failure branch",
        },
        {
            "id": "f68_failure_branch_name",
            "actual": f68["selected_branch"]["selected_branch"],
            "expected": "explicit_failure_verdict_for_construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "meaning": "F68 already fixes the explicit failure branch name",
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
            "stage": "P155",
            "lane": "current_failure_verdict_discharge_probe_for_fixed_first_contract_compliant_additive_attempt_only",
            "status": "P155_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FAILURE_VERDICT_BRANCH_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P155",
            "lane": "current_failure_verdict_discharge_probe_for_fixed_first_contract_compliant_additive_attempt_only",
            "goal": "test_whether_the_current_repo_already_exports_an_explicit_failure_verdict_discharge_for_the_fixed_first_contract_compliant_future_attempt",
            "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_FAILURE_VERDICT_DISCHARGE_FOR_THE_FIXED_FIRST_CONTRACT_COMPLIANT_FUTURE_ADDITIVE_UPSTREAM_SOURCE_CONSTRUCTION_ATTEMPT_AFTER_P155",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
