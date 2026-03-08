#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p151_first_contract_compliant_additive_upstream_source_target_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f64 = load_json(
        "fundamental_action_reconstruction/generated/f64_first_contract_compliant_additive_upstream_source_target_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "f64_target_active",
            "actual": f64["target"][
                "first_contract_compliant_additive_upstream_source_target_active"
            ],
            "expected": True,
            "meaning": "F64 already activates the first future target",
        },
        {
            "id": "f64_first_target",
            "actual": f64["target"]["first_target"],
            "expected": "S_sel_int_future_additive_upstream_target_v2",
            "meaning": "F64 already freezes one fixed first target name",
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
            "stage": "P151",
            "lane": "first_contract_compliant_additive_upstream_source_target_probe_only",
            "status": "P151_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_CONTRACT_COMPLIANT_TARGET_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P151",
            "lane": "first_contract_compliant_additive_upstream_source_target_probe_only",
            "goal": "test_whether_the_only_remaining_honest_positive_work_is_now_reduced_to_one_first_contract_compliant_future_additive_upstream_target",
            "status": "CURRENT_REPO_REDUCES_THE_ONLY_HONEST_POSITIVE_WORK_TO_ONE_FIRST_CONTRACT_COMPLIANT_FUTURE_ADDITIVE_UPSTREAM_SOURCE_TARGET_AFTER_P151",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
