#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n185_current_first_clause_obstruction_theorem_for_s_prelm_additive_candidate_v1_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f76 = load_json(
        "fundamental_action_reconstruction/generated/f76_first_additive_preobserver_source_object_construction_attempt_packet_summary.json"
    )
    n183 = load_json(
        "fundamental_action_reconstruction/generated/n183_current_first_additive_preobserver_source_object_admissibility_upgrade_target_theorem_summary.json"
    )
    n184 = load_json(
        "fundamental_action_reconstruction/generated/n184_current_first_additive_preobserver_source_object_first_clause_admissibility_target_theorem_summary.json"
    )
    p166 = load_json(
        "fundamental_action_reconstruction/generated/p166_genuinely_new_strict_core_source_object_clause_probe_for_s_prelm_additive_candidate_v1_summary.json"
    )

    checks_spec = [
        {
            "id": "f76_attempt_present",
            "actual": f76["construction_attempt"],
            "expected": "S_preLM_additive_candidate_v1",
        },
        {
            "id": "n183_upgrade_target_present",
            "actual": n183["theorem_result"]["admissibility_upgrade_target"],
            "expected": "upgrade_to_admissible_source_v1(S_preLM_additive_candidate_v1)",
        },
        {
            "id": "n184_first_clause_present",
            "actual": n184["theorem_result"]["first_clause"],
            "expected": "genuinely_new_strict_core_source_object_required",
        },
        {
            "id": "p166_clause_obstructed",
            "actual": p166["status"],
            "expected": "CURRENT_REPO_DOES_NOT_YET_SHOW_THAT_S_PRELM_ADDITIVE_CANDIDATE_V1_SATISFIES_THE_GENUINELY_NEW_STRICT_CORE_SOURCE_OBJECT_CLAUSE_AFTER_P166",
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
            "step": "N185",
            "status": "N185_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_CLAUSE_OBSTRUCTION_STATE_FOR_S_PRELM_ADDITIVE_CANDIDATE_V1",
            "scope": "current_first_clause_obstruction_for_S_preLM_additive_candidate_v1_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N185",
            "status": "N185_DISCHARGED_CURRENT_FIRST_CLAUSE_OBSTRUCTION_THEOREM_FOR_S_PRELM_ADDITIVE_CANDIDATE_V1_NO_FALSE_PASS",
            "scope": "current_first_clause_obstruction_for_S_preLM_additive_candidate_v1_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "construction_attempt": "S_preLM_additive_candidate_v1",
                "first_clause": "genuinely_new_strict_core_source_object_required",
                "full_closure_pass": False,
            },
            "hard_limits": [
                "first_clause_not_yet_satisfied",
                "later_clauses_not_yet_tested",
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
