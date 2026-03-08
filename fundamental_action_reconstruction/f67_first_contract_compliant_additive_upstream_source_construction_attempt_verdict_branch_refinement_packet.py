#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f67_first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_refinement_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f66 = load_json(
        "fundamental_action_reconstruction/generated/f66_first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_packet_summary.json"
    )
    n169 = load_json(
        "fundamental_action_reconstruction/generated/n169_current_first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f66_verdict_target_active",
            "actual": f66["verdict_target"][
                "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_active"
            ],
            "expected": True,
            "meaning": "F66 already activates the fixed verdict target",
        },
        {
            "id": "n169_verdict_target_active",
            "actual": n169["theorem_result"][
                "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_active"
            ],
            "expected": True,
            "meaning": "N169 already fixes the verdict target theorem-level",
        },
        {
            "id": "n169_verdict_target_name",
            "actual": n169["theorem_result"]["verdict_target"],
            "expected": "success_or_failure_verdict(construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2))",
            "meaning": "N169 already fixes the verdict target name",
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
            "stage": "F67",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_split_only",
            "status": "F67_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_VERDICT_BRANCH_SPLIT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F67",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_split_only",
            "goal": "freeze_the_minimal_success_failure_branch_split_inside_the_fixed_first_contract_compliant_future_attempt_verdict_target",
            "status": "F67_EXECUTED_FIRST_CONTRACT_COMPLIANT_ADDITIVE_UPSTREAM_SOURCE_CONSTRUCTION_ATTEMPT_VERDICT_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
            "checks": checks,
            "branch_split": {
                "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_split_active": True,
                "allowed_branches": [
                    "explicit_success_verdict_for_construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
                    "explicit_failure_verdict_for_construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
                ],
            },
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
