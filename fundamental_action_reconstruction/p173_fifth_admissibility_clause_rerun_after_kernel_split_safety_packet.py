#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p173_fifth_admissibility_clause_rerun_after_kernel_split_safety_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n188 = load_json(
        "fundamental_action_reconstruction/generated/n188_current_first_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    n189 = load_json(
        "fundamental_action_reconstruction/generated/n189_current_second_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    n190 = load_json(
        "fundamental_action_reconstruction/generated/n190_current_third_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    n191 = load_json(
        "fundamental_action_reconstruction/generated/n191_current_fourth_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    f85 = load_json(
        "fundamental_action_reconstruction/generated/f85_first_exported_preobserver_source_object_kernel_split_safety_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "first_clause_discharged",
            "actual": n188["theorem_result"]["discharged"],
            "expected": True,
        },
        {
            "id": "second_clause_discharged",
            "actual": n189["theorem_result"]["discharged"],
            "expected": True,
        },
        {
            "id": "third_clause_discharged",
            "actual": n190["theorem_result"]["discharged"],
            "expected": True,
        },
        {
            "id": "fourth_clause_discharged",
            "actual": n191["theorem_result"]["discharged"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe",
            "actual": f85["kernel_split_safe"],
            "expected": True,
        },
        {
            "id": "no_external_selector_import",
            "actual": f85["no_external_selector_import"],
            "expected": True,
        },
        {
            "id": "guardrail_kernel_split_safe",
            "actual": f85["guardrail_kernel_split_safe"],
            "expected": True,
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
            "stage": "P173",
            "lane": "fifth_admissibility_clause_rerun_after_kernel_split_safety_only",
            "status": "P173_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIFTH_CLAUSE_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P173",
            "lane": "fifth_admissibility_clause_rerun_after_kernel_split_safety_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_FIFTH_ADMISSIBILITY_CLAUSE_AFTER_P173",
            "exported_object": "S_preLM_strict_core_source_object_v1",
            "fifth_clause": "non_substitutive_wrt_kernel_split",
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [
                "selector_acceptance_independent",
                "future_bridge_compatible",
            ],
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
