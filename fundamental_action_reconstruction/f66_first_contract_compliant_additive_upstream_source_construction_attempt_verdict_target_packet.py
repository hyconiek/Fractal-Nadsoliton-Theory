#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f66_first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f65 = load_json(
        "fundamental_action_reconstruction/generated/f65_first_contract_compliant_additive_upstream_source_construction_attempt_packet_summary.json"
    )
    n168 = load_json(
        "fundamental_action_reconstruction/generated/n168_current_first_contract_compliant_additive_upstream_source_construction_attempt_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f65_attempt_active",
            "actual": f65["attempt"][
                "first_contract_compliant_additive_upstream_source_construction_attempt_active"
            ],
            "expected": True,
            "meaning": "F65 already activates the first contract-compliant future attempt",
        },
        {
            "id": "n168_attempt_active",
            "actual": n168["theorem_result"][
                "first_contract_compliant_additive_upstream_source_construction_attempt_active"
            ],
            "expected": True,
            "meaning": "N168 already fixes the first contract-compliant future attempt theorem-level",
        },
        {
            "id": "n168_attempt_name",
            "actual": n168["theorem_result"]["first_attempt"],
            "expected": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "meaning": "N168 already fixes the attempt name",
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
            "stage": "F66",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_only",
            "status": "F66_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_VERDICT_TARGET_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F66",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_only",
            "goal": "freeze_one_explicit_verdict_target_for_the_first_fixed_contract_compliant_future_construction_attempt",
            "status": "F66_EXECUTED_FIRST_CONTRACT_COMPLIANT_ADDITIVE_UPSTREAM_SOURCE_CONSTRUCTION_ATTEMPT_VERDICT_TARGET_PACKET_NO_FALSE_PASS",
            "checks": checks,
            "verdict_target": {
                "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_active": True,
                "verdict_target": "success_or_failure_verdict(construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2))",
            },
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
