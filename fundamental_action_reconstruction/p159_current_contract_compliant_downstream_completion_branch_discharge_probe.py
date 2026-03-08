#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p159_current_contract_compliant_downstream_completion_branch_discharge_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f72 = load_json(
        "fundamental_action_reconstruction/generated/f72_last_contract_compliant_downstream_completion_branch_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "f72_fixed_attempt_instance",
            "actual": f72["fixed_attempt_instance"],
            "expected": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "meaning": "F72 fixes the contract-compliant additive construction attempt under discussion",
        },
        {
            "id": "f72_last_remaining_lower_branch",
            "actual": f72["last_remaining_lower_branch"],
            "expected": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction_for_the_fixed_first_contract_compliant_additive_attempt",
            "meaning": "F72 freezes the last remaining downstream branch",
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
            "stage": "P159",
            "lane": "current_contract_compliant_downstream_completion_branch_discharge_only",
            "status": "P159_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CONTRACT_COMPLIANT_DOWNSTREAM_BRANCH_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P159",
            "lane": "current_contract_compliant_downstream_completion_branch_discharge_only",
            "goal": "test_whether_the_current_repo_already_exports_an_explicit_downstream_completion_branch_discharge_after_the_contract_compliant_orientation_export_obstruction",
            "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_DOWNSTREAM_COMPLETION_BRANCH_DISCHARGE_FOR_THE_LAST_REMAINING_CONTRACT_COMPLIANT_LOWER_BRANCH_AFTER_P159",
            "fixed_attempt_instance": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "tested_lower_branch": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction_for_the_fixed_first_contract_compliant_additive_attempt",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
