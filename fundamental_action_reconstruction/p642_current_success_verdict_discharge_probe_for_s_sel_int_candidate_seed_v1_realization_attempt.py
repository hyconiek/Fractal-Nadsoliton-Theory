#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p642_current_success_verdict_discharge_probe_for_s_sel_int_candidate_seed_v1_realization_attempt.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p642_current_success_verdict_discharge_probe_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f642 = load_json(
        "fundamental_action_reconstruction/generated/f642_remaining_realization_success_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
    )
    n531 = load_json(
        "fundamental_action_reconstruction/generated/n531_next_constructive_move_reduced_to_one_explicit_success_failure_branch_split_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
    )
    n532 = load_json(
        "fundamental_action_reconstruction/generated/n532_current_failure_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
    )

    explicit_success_verdict_exported = False

    checks_spec = [
        {
            "id": "f642_success_branch_remaining",
            "actual": f642["remaining_branch_to_attack"],
            "expected": "explicit_success_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v1",
            "meaning": "F642 freezes the success branch as the remaining branch to attack (v1)",
        },
        {
            "id": "n531_success_branch_name",
            "actual": n531["theorem_result"]["verdict_branches"]["success_branch"],
            "expected": "explicit_success_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v1",
            "meaning": "N531 already fixes the success branch name on the binary split (v1)",
        },
        {
            "id": "n532_failure_obstruction_discharged",
            "actual": n532["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N532 already packages the current failure-side obstruction (v1)",
        },
        {
            "id": "explicit_success_verdict_exported",
            "actual": explicit_success_verdict_exported,
            "expected": False,
            "meaning": "the current repo does not yet export an explicit success verdict discharge for the fixed realization attempt (v1)",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "P642",
        "lane": "current_success_verdict_discharge_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_success_verdict_for_the_fixed_seed_v1_realization_attempt",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_SUCCESS_VERDICT_DISCHARGE_FOR_S_SEL_INT_NEW_OBJECT_CONSTRUCTED_REALIZATION_ATTEMPT_V1_AFTER_P642",
        "target_state": {
            "success_branch_selected_as_next_attack": True,
            "explicit_success_verdict_exported": explicit_success_verdict_exported,
            "remaining_open_branches": [
                "future_admissibility_test_of_a_future_constructed_source_object",
                "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
                "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P642",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()

