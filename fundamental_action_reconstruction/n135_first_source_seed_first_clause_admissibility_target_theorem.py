#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n135_first_source_seed_first_clause_admissibility_target_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p124 = load_json(
        "fundamental_action_reconstruction/generated/p124_first_source_seed_genuinely_new_strict_core_object_clause_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p124_first_clause_reduced",
            "actual": p124["target_state"]["next_clause_by_clause_move_reduced_to_first_clause"],
            "expected": True,
            "meaning": "P124 already proves the next clause-by-clause move is reduced to the first clause",
        },
        {
            "id": "first_clause_name",
            "actual": p124["target_state"]["first_clause_target"]["first_clause_name"],
            "expected": "genuinely_new_strict_core_source_object_required",
            "meaning": "the first clause is explicitly fixed",
        },
        {
            "id": "candidate_seed_name",
            "actual": p124["target_state"]["first_clause_target"]["candidate_seed_name"],
            "expected": "S_sel_int_candidate_seed_v0",
            "meaning": "the first clause is anchored on the one candidate seed instance",
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
            "step": "N135",
            "status": "N135_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_CLAUSE_TARGET_STATE",
            "scope": "first_source_seed_first_clause_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N135",
            "status": "N135_DISCHARGED_FIRST_SOURCE_SEED_FIRST_CLAUSE_ADMISSIBILITY_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "first_source_seed_first_clause_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_clause_by_clause_move_reduced_to_first_clause": True,
                "first_clause_target": p124["target_state"]["first_clause_target"],
                "later_clauses_left_open": p124["target_state"]["later_clauses_left_open"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_test_of_whether_S_sel_int_candidate_seed_v0_counts_as_a_genuinely_new_strict_core_source_object",
                "future_clause_by_clause_tests_for_the_remaining_F34_admissibility_requirements",
                "future_derivation_of_admissible_E_orient_from_an_upgraded_source_seed",
                "future_completion_of_B_sel_R_sel_O_sel_after_seed_upgrade",
            ],
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
