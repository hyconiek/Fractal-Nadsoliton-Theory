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
    / "p643_current_admissibility_branch_discharge_probe_for_future_constructed_source_object_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p643_current_admissibility_branch_discharge_probe_for_future_constructed_source_object_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f643 = load_json(
        "fundamental_action_reconstruction/generated/f643_first_post_verdict_admissibility_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
    )
    n532 = load_json(
        "fundamental_action_reconstruction/generated/n532_current_failure_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
    )
    n533 = load_json(
        "fundamental_action_reconstruction/generated/n533_current_success_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
    )

    explicit_admissibility_branch_discharge_exported = False

    checks_spec = [
        {
            "id": "f643_admissibility_branch_first",
            "actual": f643["first_lower_branch_to_attack"],
            "expected": "future_admissibility_test_of_a_future_constructed_source_object",
            "meaning": "F643 freezes admissibility as the first remaining lower branch (v1)",
        },
        {
            "id": "n532_failure_obstruction_discharged",
            "actual": n532["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N532 packages the failure-side obstruction (v1)",
        },
        {
            "id": "n533_success_obstruction_discharged",
            "actual": n533["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N533 packages the success-side obstruction (v1)",
        },
        {
            "id": "explicit_admissibility_branch_discharge_exported",
            "actual": explicit_admissibility_branch_discharge_exported,
            "expected": False,
            "meaning": "the current repo does not export an explicit admissibility-branch discharge for seed-v1 at this scope",
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
        "stage": "P643",
        "lane": "current_admissibility_branch_discharge_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_admissibility_branch_discharge_for_seed_v1_post_verdict_admissibility_branch",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_ADMISSIBILITY_BRANCH_DISCHARGE_FOR_SEED_V1_AFTER_P643",
        "target_state": {
            "admissibility_branch_selected_as_next_attack": True,
            "explicit_admissibility_branch_discharge_exported": explicit_admissibility_branch_discharge_exported,
            "remaining_open_branches": [
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
        "stage": "P643",
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

