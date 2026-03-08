#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p154_first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f67 = load_json(
        "fundamental_action_reconstruction/generated/f67_first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_refinement_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "f67_branch_split_active",
            "actual": f67["branch_split"][
                "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_split_active"
            ],
            "expected": True,
            "meaning": "F67 already activates the explicit branch split",
        },
        {
            "id": "f67_allowed_branch_count",
            "actual": len(f67["branch_split"]["allowed_branches"]),
            "expected": 2,
            "meaning": "F67 already reduces the verdict target to exactly two branches",
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
            "stage": "P154",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_probe_only",
            "status": "P154_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_VERDICT_BRANCH_SPLIT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P154",
            "lane": "first_contract_compliant_additive_upstream_source_construction_attempt_verdict_branch_probe_only",
            "goal": "test_whether_the_fixed_first_contract_compliant_future_attempt_verdict_target_is_now_reduced_to_one_explicit_success_failure_branch_split",
            "status": "CURRENT_REPO_REDUCES_THE_FIXED_FIRST_CONTRACT_COMPLIANT_FUTURE_ADDITIVE_UPSTREAM_SOURCE_CONSTRUCTION_ATTEMPT_VERDICT_TARGET_TO_ONE_EXPLICIT_SUCCESS_FAILURE_BRANCH_SPLIT_AFTER_P154",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
