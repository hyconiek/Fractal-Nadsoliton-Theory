#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n184_current_first_additive_preobserver_source_object_first_clause_admissibility_target_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n183 = load_json(
        "fundamental_action_reconstruction/generated/n183_current_first_additive_preobserver_source_object_admissibility_upgrade_target_theorem_summary.json"
    )
    f78 = load_json(
        "fundamental_action_reconstruction/generated/f78_first_additive_preobserver_source_object_first_clause_refinement_packet_summary.json"
    )
    p165 = load_json(
        "fundamental_action_reconstruction/generated/p165_first_additive_preobserver_source_object_genuinely_new_strict_core_object_clause_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "n183_upgrade_target_present",
            "actual": n183["theorem_result"]["admissibility_upgrade_target"],
            "expected": "upgrade_to_admissible_source_v1(S_preLM_additive_candidate_v1)",
        },
        {
            "id": "f78_first_clause_present",
            "actual": f78["first_clause"],
            "expected": "genuinely_new_strict_core_source_object_required",
        },
        {
            "id": "p165_clause_probe_positive",
            "actual": p165["first_clause"],
            "expected": "genuinely_new_strict_core_source_object_required",
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
            "step": "N184",
            "status": "N184_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_ADDITIVE_PREOBSERVER_FIRST_CLAUSE_STATE",
            "scope": "current_first_additive_preobserver_source_object_first_clause_admissibility_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N184",
            "status": "N184_DISCHARGED_CURRENT_FIRST_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_FIRST_CLAUSE_ADMISSIBILITY_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "current_first_additive_preobserver_source_object_first_clause_admissibility_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "construction_attempt": "S_preLM_additive_candidate_v1",
                "admissibility_upgrade_target": "upgrade_to_admissible_source_v1(S_preLM_additive_candidate_v1)",
                "first_clause": "genuinely_new_strict_core_source_object_required",
                "full_closure_pass": False,
            },
            "hard_limits": [
                "first_clause_not_yet_satisfied",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
