#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p175_seventh_admissibility_clause_rerun_after_future_bridge_compatibility_packet_summary.json"


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
    n192 = load_json(
        "fundamental_action_reconstruction/generated/n192_current_fifth_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    n193 = load_json(
        "fundamental_action_reconstruction/generated/n193_current_sixth_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    f87 = load_json(
        "fundamental_action_reconstruction/generated/f87_first_exported_preobserver_source_object_future_bridge_compatibility_packet_summary.json"
    )

    checks_spec = [
        {"id": "first_clause_discharged", "actual": n188["theorem_result"]["discharged"], "expected": True},
        {"id": "second_clause_discharged", "actual": n189["theorem_result"]["discharged"], "expected": True},
        {"id": "third_clause_discharged", "actual": n190["theorem_result"]["discharged"], "expected": True},
        {"id": "fourth_clause_discharged", "actual": n191["theorem_result"]["discharged"], "expected": True},
        {"id": "fifth_clause_discharged", "actual": n192["theorem_result"]["discharged"], "expected": True},
        {"id": "sixth_clause_discharged", "actual": n193["theorem_result"]["discharged"], "expected": True},
        {"id": "source_object_first", "actual": f87["source_object_first"], "expected": True},
        {"id": "upstream_of_observer", "actual": f87["upstream_of_observer"], "expected": True},
        {"id": "second_clause_bridge_ready", "actual": f87["second_clause_discharged"], "expected": True},
        {"id": "third_clause_bridge_ready", "actual": f87["third_clause_discharged"], "expected": True},
        {"id": "fifth_clause_bridge_ready", "actual": f87["fifth_clause_discharged"], "expected": True},
        {"id": "sixth_clause_bridge_ready", "actual": f87["sixth_clause_discharged"], "expected": True},
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append({"id": item["id"], "actual": item["actual"], "expected": item["expected"], "pass": ok})
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "stage": "P175",
            "lane": "seventh_admissibility_clause_rerun_after_future_bridge_compatibility_only",
            "status": "P175_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SEVENTH_CLAUSE_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P175",
            "lane": "seventh_admissibility_clause_rerun_after_future_bridge_compatibility_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_SEVENTH_ADMISSIBILITY_CLAUSE_AFTER_P175",
            "exported_object": "S_preLM_strict_core_source_object_v1",
            "seventh_clause": "future_bridge_compatible",
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [],
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
