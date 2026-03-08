#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p171_third_admissibility_clause_rerun_after_source_seed_only_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n188 = load_json(
        "fundamental_action_reconstruction/generated/n188_current_first_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    n189 = load_json(
        "fundamental_action_reconstruction/generated/n189_current_second_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    f83 = load_json(
        "fundamental_action_reconstruction/generated/f83_first_exported_preobserver_source_object_source_seed_only_packet_summary.json"
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
            "id": "source_seed_only",
            "actual": f83["source_seed_only"],
            "expected": True,
        },
        {
            "id": "E_orient_not_exported",
            "actual": f83["E_orient_exported"],
            "expected": False,
        },
        {
            "id": "B_sel_not_exported",
            "actual": f83["B_sel_exported"],
            "expected": False,
        },
        {
            "id": "R_sel_not_exported",
            "actual": f83["R_sel_exported"],
            "expected": False,
        },
        {
            "id": "O_sel_not_exported",
            "actual": f83["O_sel_exported"],
            "expected": False,
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
            "stage": "P171",
            "lane": "third_admissibility_clause_rerun_after_source_seed_only_only",
            "status": "P171_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_THIRD_CLAUSE_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P171",
            "lane": "third_admissibility_clause_rerun_after_source_seed_only_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_THIRD_ADMISSIBILITY_CLAUSE_AFTER_P171",
            "exported_object": "S_preLM_strict_core_source_object_v1",
            "third_clause": "source_seed_only",
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [
                "non_substitutive",
                "selector_acceptance_independent",
                "future_bridge_compatible",
            ],
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
