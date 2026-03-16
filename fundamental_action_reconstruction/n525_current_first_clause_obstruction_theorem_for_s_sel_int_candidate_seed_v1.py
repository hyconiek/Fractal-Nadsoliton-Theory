#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n525_current_first_clause_obstruction_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p634 = load_json(
        "fundamental_action_reconstruction/generated/p634_genuinely_new_strict_core_source_object_clause_probe_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    checks_spec = [
        {
            "id": "p634_first_clause_currently_unsatisfied",
            "actual": p634["clause_test_result"]["currently_satisfied"],
            "expected": False,
            "meaning": "P634 shows the first clause is not yet satisfied for seed v1",
        },
        {
            "id": "first_clause_name",
            "actual": p634["clause_test_result"]["first_clause_name"],
            "expected": "genuinely_new_strict_core_source_object_required",
            "meaning": "the theorem stays scoped to the first admissibility clause",
        },
        {
            "id": "candidate_seed_name",
            "actual": p634["clause_test_result"]["candidate_seed_name"],
            "expected": "S_sel_int_candidate_seed_v1",
            "meaning": "the theorem stays scoped to the strict sigma-int upgraded seed candidate instance",
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
            "step": "N525",
            "status": "N525_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_CLAUSE_OBSTRUCTION_STATE",
            "scope": "first_source_seed_first_clause_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
    else:
        summary = {
            "step": "N525",
            "status": "N525_DISCHARGED_CURRENT_FIRST_CLAUSE_OBSTRUCTION_THEOREM_FOR_S_SEL_INT_CANDIDATE_SEED_V1_NO_FALSE_PASS",
            "scope": "first_source_seed_first_clause_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "candidate_seed_name": "S_sel_int_candidate_seed_v1",
                "first_clause_name": "genuinely_new_strict_core_source_object_required",
                "first_clause_currently_satisfied": False,
                "strict_sigma_int_slot_upgraded": bool(
                    p634["clause_test_result"].get("strict_sigma_int_slot_upgraded")
                ),
                "current_reason": p634["clause_test_result"]["current_reason"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_attempt_to_construct_a_genuinely_new_strict_core_source_object_beyond_the_upgraded_seed_v1",
                "future_clause_by_clause_tests_for_the_remaining_F34_admissibility_requirements",
                "future_derivation_of_admissible_E_orient_from_a_future_upgraded_source_object",
                "future_completion_of_B_sel_R_sel_O_sel_after_source_object_construction",
            ],
            "hard_limits": [
                "first_clause_currently_unsatisfied",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

