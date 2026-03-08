#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p165_first_additive_preobserver_source_object_genuinely_new_strict_core_object_clause_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f77 = load_json(
        "fundamental_action_reconstruction/generated/f77_first_additive_preobserver_source_object_admissibility_upgrade_target_packet_summary.json"
    )
    p164 = load_json(
        "fundamental_action_reconstruction/generated/p164_current_additive_preobserver_source_object_admissibility_upgrade_target_probe_summary.json"
    )
    n183 = load_json(
        "fundamental_action_reconstruction/generated/n183_current_first_additive_preobserver_source_object_admissibility_upgrade_target_theorem_summary.json"
    )
    f78 = load_json(
        "fundamental_action_reconstruction/generated/f78_first_additive_preobserver_source_object_first_clause_refinement_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "upgrade_target_present",
            "actual": f77["admissibility_upgrade_target"],
            "expected": "upgrade_to_admissible_source_v1(S_preLM_additive_candidate_v1)",
        },
        {
            "id": "upgrade_target_probe_positive",
            "actual": p164["admissibility_upgrade_target"],
            "expected": "upgrade_to_admissible_source_v1(S_preLM_additive_candidate_v1)",
        },
        {
            "id": "upgrade_target_theorem_future_only",
            "actual": [
                n183["theorem_result"]["future_only"],
                n183["theorem_result"]["upstream_of_observer"],
            ],
            "expected": [True, True],
        },
        {
            "id": "first_clause_packet_present",
            "actual": f78["first_clause"],
            "expected": "genuinely_new_strict_core_source_object_required",
        },
        {
            "id": "contract_source_preserved",
            "actual": f78["contract_source"],
            "expected": "F34",
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
            "stage": "P165",
            "lane": "future_additive_preobserver_source_object_first_clause_probe_only",
            "status": "P165_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_CLAUSE_TARGET_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P165",
            "lane": "future_additive_preobserver_source_object_first_clause_probe_only",
            "goal": "test_whether_the_first_clause_by_clause_admissibility_question_is_now_reduced_to_the_genuinely_new_strict_core_source_object_requirement",
            "status": "CURRENT_REPO_REDUCES_THE_FIRST_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_CLAUSE_BY_CLAUSE_ADMISSIBILITY_TEST_TO_THE_GENUINELY_NEW_STRICT_CORE_SOURCE_OBJECT_REQUIREMENT_AFTER_P165",
            "construction_attempt": "S_preLM_additive_candidate_v1",
            "first_clause": "genuinely_new_strict_core_source_object_required",
            "checks": checks,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
