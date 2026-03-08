#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p152_first_contract_compliant_additive_upstream_source_construction_attempt_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f65 = load_json(
        "fundamental_action_reconstruction/generated/f65_first_contract_compliant_additive_upstream_source_construction_attempt_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "f65_attempt_active",
            "actual": f65["attempt"][
                "first_contract_compliant_additive_upstream_source_construction_attempt_active"
            ],
            "expected": True,
            "meaning": "F65 already activates the first future construction attempt",
        },
        {
            "id": "f65_attempt_name",
            "actual": f65["attempt"]["attempt"],
            "expected": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "meaning": "F65 already fixes one explicit first attempt name",
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
            "stage": "P152",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_probe_only",
            "status": "P152_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_CONTRACT_COMPLIANT_CONSTRUCTION_ATTEMPT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P152",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_probe_only",
            "goal": "test_whether_the_only_remaining_honest_positive_work_is_now_reduced_to_one_first_contract_compliant_future_additive_upstream_construction_attempt",
            "status": "CURRENT_REPO_REDUCES_THE_ONLY_HONEST_POSITIVE_WORK_TO_ONE_FIRST_CONTRACT_COMPLIANT_FUTURE_ADDITIVE_UPSTREAM_SOURCE_CONSTRUCTION_ATTEMPT_AFTER_P152",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
