#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p168_current_additive_preobserver_source_object_first_clause_remainder_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f80 = load_json(
        "fundamental_action_reconstruction/generated/f80_first_additive_preobserver_source_object_first_clause_remainder_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "nonreduction_witness_present",
            "actual": f80["nonreduction_witness_present"],
            "expected": True,
        },
        {
            "id": "simple_packaging_reduction_removed",
            "actual": f80["simple_packaging_reduction_removed"],
            "expected": True,
        },
        {
            "id": "remaining_first_clause_package",
            "actual": f80["remaining_first_clause_package"],
            "expected": "realized_constructed_source_object_export_package",
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
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "stage": "P168",
            "lane": "current_additive_preobserver_source_object_first_clause_remainder_probe_only",
            "status": "P168_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_CLAUSE_REMAINDER_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P168",
            "lane": "current_additive_preobserver_source_object_first_clause_remainder_probe_only",
            "status": "CURRENT_REPO_REDUCES_THE_FIRST_CLAUSE_OBSTRUCTION_FOR_S_PRELM_ADDITIVE_CANDIDATE_V1_TO_ONE_REALIZED_CONSTRUCTED_SOURCE_OBJECT_EXPORT_PACKAGE_AFTER_P168",
            "construction_attempt": "S_preLM_additive_candidate_v1",
            "first_clause": "genuinely_new_strict_core_source_object_required",
            "remaining_first_clause_package": "realized_constructed_source_object_export_package",
            "remaining_first_clause_obstructions": f80["remaining_first_clause_obstructions"],
            "checks": checks,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
