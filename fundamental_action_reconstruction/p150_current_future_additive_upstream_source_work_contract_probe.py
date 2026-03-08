#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p150_current_future_additive_upstream_source_work_contract_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f63 = load_json(
        "fundamental_action_reconstruction/generated/f63_current_future_additive_upstream_source_work_contract_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "f63_contract_active",
            "actual": f63["contract"][
                "current_future_additive_upstream_source_work_contract_active"
            ],
            "expected": True,
            "meaning": "F63 already activates the additive upstream source-work contract",
        },
        {
            "id": "f63_kernel_split_safe",
            "actual": f63["contract"]["admitted_contract"]["kernel_split_safe"],
            "expected": True,
            "meaning": "F63 keeps the remaining positive work kernel-split safe",
        },
        {
            "id": "f63_upstream_of_observer",
            "actual": f63["contract"]["admitted_contract"]["upstream_of_observer"],
            "expected": True,
            "meaning": "F63 keeps the remaining positive work upstream of observer",
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
            "stage": "P150",
            "lane": "current_future_additive_upstream_source_work_contract_probe_only",
            "status": "P150_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADDITIVE_UPSTREAM_WORK_CONTRACT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P150",
            "lane": "current_future_additive_upstream_source_work_contract_probe_only",
            "goal": "test_whether_the_only_honest_positive_work_left_after_handoff_is_reduced_to_one_explicit_additive_upstream_contract",
            "status": "CURRENT_REPO_REDUCES_THE_ONLY_HONEST_POSITIVE_WORK_TO_ONE_EXPLICIT_FUTURE_ADDITIVE_UPSTREAM_SOURCE_WORK_CONTRACT_AFTER_P150",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
