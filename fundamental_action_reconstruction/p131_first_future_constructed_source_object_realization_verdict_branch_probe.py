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
    / "p131_first_future_constructed_source_object_realization_verdict_branch_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p131_first_future_constructed_source_object_realization_verdict_branch_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n141 = load_json(
        "fundamental_action_reconstruction/generated/n141_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_verdict_target_theorem_summary.json"
    )
    f44 = load_json(
        "fundamental_action_reconstruction/generated/f44_first_future_constructed_source_object_realization_verdict_branch_refinement_packet_summary.json"
    )

    reduced_to_binary_branch_split = (
        n141["theorem_result"]["next_constructive_move_reduced_to_one_first_future_verdict_target"]
        is True
        and f44["verdict_branches"]["success_branch"]
        == "explicit_success_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0"
        and f44["verdict_branches"]["failure_branch"]
        == "explicit_failure_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0"
    )

    checks_spec = [
        {
            "id": "n141_verdict_target_reduced",
            "actual": n141["theorem_result"][
                "next_constructive_move_reduced_to_one_first_future_verdict_target"
            ],
            "expected": True,
            "meaning": "N141 already reduces the next constructive move to one verdict target",
        },
        {
            "id": "success_branch_name",
            "actual": f44["verdict_branches"]["success_branch"],
            "expected": "explicit_success_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0",
            "meaning": "F44 exports one explicit success branch",
        },
        {
            "id": "failure_branch_name",
            "actual": f44["verdict_branches"]["failure_branch"],
            "expected": "explicit_failure_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0",
            "meaning": "F44 exports one explicit failure branch",
        },
        {
            "id": "reduced_to_binary_branch_split",
            "actual": reduced_to_binary_branch_split,
            "expected": True,
            "meaning": "the next constructive move is reduced to one explicit binary verdict split",
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
        "stage": "P131",
        "lane": "first_future_constructed_source_object_realization_verdict_branch_refinement_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_explicit_binary_success_failure_branch_split",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_EXPLICIT_SUCCESS_FAILURE_BRANCH_SPLIT_AFTER_P131",
        "target_state": {
            "next_constructive_move_reduced_to_binary_branch_split": reduced_to_binary_branch_split,
            "verdict_target": f44["verdict_target"],
            "verdict_branches": f44["verdict_branches"],
            "later_open_branches": [
                "future_success_verdict_discharge_for_S_sel_int_new_object_constructed_realization_attempt_v0",
                "future_failure_verdict_discharge_for_S_sel_int_new_object_constructed_realization_attempt_v0",
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
        "stage": "P131",
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
