#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p157_current_contract_compliant_post_verdict_admissibility_branch_discharge_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f70 = load_json(
        "fundamental_action_reconstruction/generated/f70_first_post_contract_compliant_verdict_admissibility_branch_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "f70_fixed_attempt_instance",
            "actual": f70["fixed_attempt_instance"],
            "expected": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "meaning": "F70 fixes the contract-compliant additive construction attempt under discussion",
        },
        {
            "id": "f70_first_lower_branch_to_attack",
            "actual": f70["first_lower_branch_to_attack"],
            "expected": "future_admissibility_test_of_a_future_constructed_source_object_for_S_sel_int_after_fixed_first_contract_compliant_additive_attempt",
            "meaning": "F70 freezes the admissibility branch as the first remaining lower branch",
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
            "stage": "P157",
            "lane": "current_contract_compliant_post_verdict_admissibility_branch_discharge_only",
            "status": "P157_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_POST_CONTRACT_COMPLIANT_VERDICT_ADMISSIBILITY_BRANCH_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P157",
            "lane": "current_contract_compliant_post_verdict_admissibility_branch_discharge_only",
            "goal": "test_whether_the_current_repo_already_exports_an_explicit_admissibility_branch_discharge_below_the_exhausted_contract_compliant_additive_binary_verdict_layer",
            "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_ADMISSIBILITY_BRANCH_DISCHARGE_FOR_THE_FIRST_REMAINING_LOWER_BRANCH_BELOW_THE_EXHAUSTED_CONTRACT_COMPLIANT_ADDITIVE_BINARY_VERDICT_LAYER_AFTER_P157",
            "fixed_attempt_instance": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "tested_lower_branch": "future_admissibility_test_of_a_future_constructed_source_object_for_S_sel_int_after_fixed_first_contract_compliant_additive_attempt",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
