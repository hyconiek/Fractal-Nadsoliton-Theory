#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p153_first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f66 = load_json(
        "fundamental_action_reconstruction/generated/f66_first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "f66_verdict_target_active",
            "actual": f66["verdict_target"][
                "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_active"
            ],
            "expected": True,
            "meaning": "F66 already activates the verdict target",
        },
        {
            "id": "f66_verdict_target_name",
            "actual": f66["verdict_target"]["verdict_target"],
            "expected": "success_or_failure_verdict(construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2))",
            "meaning": "F66 already fixes one explicit verdict target name",
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
            "stage": "P153",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_probe_only",
            "status": "P153_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_VERDICT_TARGET_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P153",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_target_probe_only",
            "goal": "test_whether_the_only_remaining_honest_positive_work_is_now_reduced_to_one_explicit_future_success_or_failure_verdict_target_for_the_fixed_first_attempt",
            "status": "CURRENT_REPO_REDUCES_THE_ONLY_HONEST_POSITIVE_WORK_TO_ONE_EXPLICIT_FUTURE_SUCCESS_OR_FAILURE_VERDICT_TARGET_FOR_THE_FIXED_FIRST_CONTRACT_COMPLIANT_ADDITIVE_UPSTREAM_SOURCE_CONSTRUCTION_ATTEMPT_AFTER_P153",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
