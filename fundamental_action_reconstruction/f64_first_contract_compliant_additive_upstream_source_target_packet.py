#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f64_first_contract_compliant_additive_upstream_source_target_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f63 = load_json(
        "fundamental_action_reconstruction/generated/f63_current_future_additive_upstream_source_work_contract_packet_summary.json"
    )
    n166 = load_json(
        "fundamental_action_reconstruction/generated/n166_current_only_honest_positive_work_contract_theorem_summary.json"
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
            "id": "n166_genuinely_additive",
            "actual": n166["theorem_result"][
                "remaining_positive_work_must_be_genuinely_additive"
            ],
            "expected": True,
            "meaning": "N166 already restricts remaining work to genuinely additive moves",
        },
        {
            "id": "n166_upstream_of_observer",
            "actual": n166["theorem_result"][
                "remaining_positive_work_must_stay_upstream_of_observer"
            ],
            "expected": True,
            "meaning": "N166 already keeps remaining work upstream of observer",
        },
        {
            "id": "n166_kernel_split_safe",
            "actual": n166["theorem_result"][
                "remaining_positive_work_must_be_kernel_split_safe"
            ],
            "expected": True,
            "meaning": "N166 already keeps remaining work kernel-split safe",
        },
        {
            "id": "n166_source_object_first",
            "actual": n166["theorem_result"][
                "remaining_positive_work_must_be_source_object_first"
            ],
            "expected": True,
            "meaning": "N166 already keeps remaining work source-object-first",
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
            "stage": "F64",
            "lane": "first_contract_compliant_additive_upstream_source_target_only",
            "status": "F64_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CONTRACT_COMPLIANT_TARGET_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F64",
            "lane": "first_contract_compliant_additive_upstream_source_target_only",
            "goal": "freeze_one_first_explicit_future_target_under_the_only_remaining_honest_positive_work_contract",
            "status": "F64_EXECUTED_FIRST_CONTRACT_COMPLIANT_ADDITIVE_UPSTREAM_SOURCE_TARGET_PACKET_NO_FALSE_PASS",
            "checks": checks,
            "target": {
                "first_contract_compliant_additive_upstream_source_target_active": True,
                "first_target": "S_sel_int_future_additive_upstream_target_v2",
            },
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
