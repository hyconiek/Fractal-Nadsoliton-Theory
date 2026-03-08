#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f72_last_contract_compliant_downstream_completion_branch_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n173 = load_json(
        "fundamental_action_reconstruction/generated/n173_current_contract_compliant_post_verdict_admissibility_branch_obstruction_theorem_summary.json"
    )
    n174 = load_json(
        "fundamental_action_reconstruction/generated/n174_current_contract_compliant_orientation_export_branch_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n173_admissibility_obstructed",
            "actual": n173["theorem_result"]["current_contract_compliant_post_verdict_admissibility_branch_obstructed"],
            "expected": True,
            "meaning": "N173 already obstructs the contract-compliant admissibility branch",
        },
        {
            "id": "n174_orientation_export_obstructed",
            "actual": n174["theorem_result"]["current_contract_compliant_orientation_export_branch_obstructed"],
            "expected": True,
            "meaning": "N174 already obstructs the contract-compliant orientation-export branch",
        },
        {
            "id": "fixed_attempt_match",
            "actual": n174["theorem_result"]["fixed_attempt_instance"],
            "expected": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "meaning": "N174 already fixes the same contract-compliant additive attempt",
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
            "stage": "F72",
            "lane": "last_contract_compliant_downstream_completion_branch_only",
            "status": "F72_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CONTRACT_COMPLIANT_DOWNSTREAM_BRANCH_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F72",
            "lane": "last_contract_compliant_downstream_completion_branch_only",
            "goal": "freeze_the_contract_compliant_downstream_completion_branch_as_the_only_remaining_lower_branch",
            "status": "F72_EXECUTED_LAST_CONTRACT_COMPLIANT_DOWNSTREAM_COMPLETION_BRANCH_PACKET_NO_FALSE_PASS",
            "fixed_attempt_instance": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "last_remaining_lower_branch": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction_for_the_fixed_first_contract_compliant_additive_attempt",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
