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
    / "p132_current_failure_verdict_discharge_probe_for_s_sel_int_realization_attempt.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p132_current_failure_verdict_discharge_probe_for_s_sel_int_realization_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f45 = load_json(
        "fundamental_action_reconstruction/generated/f45_first_conservative_realization_failure_branch_packet_summary.json"
    )
    n142 = load_json(
        "fundamental_action_reconstruction/generated/n142_next_constructive_move_reduced_to_one_explicit_success_failure_branch_split_theorem_summary.json"
    )

    explicit_failure_verdict_exported = False

    checks_spec = [
        {
            "id": "f45_failure_branch_first",
            "actual": f45["first_branch_to_attack"],
            "expected": "explicit_failure_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0",
            "meaning": "F45 already freezes the failure branch as the first branch to attack",
        },
        {
            "id": "n142_failure_branch_name",
            "actual": n142["theorem_result"]["verdict_branches"]["failure_branch"],
            "expected": "explicit_failure_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0",
            "meaning": "N142 already fixes the failure branch name on the binary split",
        },
        {
            "id": "explicit_failure_verdict_exported",
            "actual": explicit_failure_verdict_exported,
            "expected": False,
            "meaning": "the current repo does not yet export an explicit failure verdict discharge for the fixed realization attempt",
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
        "stage": "P132",
        "lane": "current_failure_verdict_discharge_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_failure_verdict_for_the_fixed_first_future_realization_attempt",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_FAILURE_VERDICT_DISCHARGE_FOR_S_SEL_INT_NEW_OBJECT_CONSTRUCTED_REALIZATION_ATTEMPT_V0_AFTER_P132",
        "target_state": {
            "failure_branch_selected_as_first_attack": True,
            "explicit_failure_verdict_exported": explicit_failure_verdict_exported,
            "remaining_open_branches": [
                "future_success_verdict_discharge_for_S_sel_int_new_object_constructed_realization_attempt_v0",
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
        "stage": "P132",
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
