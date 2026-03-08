#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f70_first_post_contract_compliant_verdict_admissibility_branch_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n171 = load_json(
        "fundamental_action_reconstruction/generated/n171_current_failure_branch_obstruction_theorem_for_fixed_first_contract_compliant_additive_attempt_summary.json"
    )
    n172 = load_json(
        "fundamental_action_reconstruction/generated/n172_current_success_branch_obstruction_theorem_for_fixed_first_contract_compliant_additive_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "n171_failure_branch_obstructed",
            "actual": n171["theorem_result"]["failure_branch_current_state_obstructed"],
            "expected": True,
            "meaning": "N171 already closes the failure branch negatively on the current repo state",
        },
        {
            "id": "n172_success_branch_obstructed",
            "actual": n172["theorem_result"]["success_branch_current_state_obstructed"],
            "expected": True,
            "meaning": "N172 already closes the success branch negatively on the current repo state",
        },
        {
            "id": "n172_binary_layer_negative",
            "actual": n172["theorem_result"]["binary_verdict_layer_fully_negative_on_current_repo_state"],
            "expected": True,
            "meaning": "N172 already exhausts the binary verdict layer for the fixed attempt",
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
            "stage": "F70",
            "lane": "post_contract_compliant_verdict_branch_order_only",
            "status": "F70_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_POST_CONTRACT_COMPLIANT_VERDICT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F70",
            "lane": "post_contract_compliant_verdict_admissibility_first",
            "goal": "freeze_the_admissibility_branch_as_the_first_remaining_lower_branch_below_the_exhausted_contract_compliant_additive_verdict_layer",
            "status": "F70_EXECUTED_FIRST_POST_CONTRACT_COMPLIANT_VERDICT_ADMISSIBILITY_BRANCH_PACKET_NO_FALSE_PASS",
            "fixed_attempt_instance": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "first_lower_branch_to_attack": "future_admissibility_test_of_a_future_constructed_source_object_for_S_sel_int_after_fixed_first_contract_compliant_additive_attempt",
            "basis": [
                "binary_verdict_layer_already_closed_by_N171_and_N172",
                "E_orient_still_presupposes_an_admissible_source_object",
                "downstream_B_sel_R_sel_O_sel_still_presupposes_source_and_export_stages",
                "preferred_order_remains_nadsoliton_to_light_to_matter_to_emergent_observer",
            ],
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
