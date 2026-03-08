#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f65_first_contract_compliant_additive_upstream_source_construction_attempt_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f63 = load_json(
        "fundamental_action_reconstruction/generated/f63_current_future_additive_upstream_source_work_contract_packet_summary.json"
    )
    n167 = load_json(
        "fundamental_action_reconstruction/generated/n167_current_first_contract_compliant_additive_upstream_source_target_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f63_contract_active",
            "actual": f63["contract"][
                "current_future_additive_upstream_source_work_contract_active"
            ],
            "expected": True,
            "meaning": "F63 already activates the additive upstream contract",
        },
        {
            "id": "n167_target_active",
            "actual": n167["theorem_result"][
                "first_contract_compliant_additive_upstream_source_target_active"
            ],
            "expected": True,
            "meaning": "N167 already activates one first contract-compliant future target",
        },
        {
            "id": "n167_target_name",
            "actual": n167["theorem_result"]["first_target"],
            "expected": "S_sel_int_future_additive_upstream_target_v2",
            "meaning": "N167 already fixes the first target name",
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
            "stage": "F65",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_only",
            "status": "F65_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_CONTRACT_COMPLIANT_CONSTRUCTION_ATTEMPT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F65",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_only",
            "goal": "freeze_one_first_future_construction_attempt_under_the_explicit_additive_upstream_contract",
            "status": "F65_EXECUTED_FIRST_CONTRACT_COMPLIANT_ADDITIVE_UPSTREAM_SOURCE_CONSTRUCTION_ATTEMPT_PACKET_NO_FALSE_PASS",
            "checks": checks,
            "attempt": {
                "first_contract_compliant_additive_upstream_source_construction_attempt_active": True,
                "attempt": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            },
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
