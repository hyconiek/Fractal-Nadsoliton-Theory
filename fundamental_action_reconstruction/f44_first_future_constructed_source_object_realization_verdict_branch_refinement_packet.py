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
    / "f44_first_future_constructed_source_object_realization_verdict_branch_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f44_first_future_constructed_source_object_realization_verdict_branch_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n141 = load_json(
        "fundamental_action_reconstruction/generated/n141_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_verdict_target_theorem_summary.json"
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
            "id": "n141_verdict_shape",
            "actual": n141["theorem_result"]["first_future_realization_verdict_target"][
                "verdict_shape"
            ],
            "expected": "success_or_failure_verdict",
            "meaning": "the fixed verdict target is already binary at this scope",
        },
        {
            "id": "n141_target_name",
            "actual": n141["theorem_result"]["first_future_realization_verdict_target"][
                "target_name"
            ],
            "expected": "S_sel_int_new_object_constructed_realization_verdict_target_v0",
            "meaning": "the verdict target name is fixed",
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

    verdict_branches = {
        "success_branch": "explicit_success_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0",
        "failure_branch": "explicit_failure_verdict_for_S_sel_int_new_object_constructed_realization_attempt_v0",
    }

    artifact = {
        "stage": "F44",
        "lane": "first_future_constructed_source_object_realization_verdict_branch_refinement_only",
        "goal": "freeze_the_two_explicit_verdict_branches_on_the_fixed_realization_verdict_target",
        "status": "F44_EXECUTED_FIRST_FUTURE_CONSTRUCTED_SOURCE_OBJECT_REALIZATION_VERDICT_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
        "verdict_target": n141["theorem_result"]["first_future_realization_verdict_target"],
        "verdict_branches": verdict_branches,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F44",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "verdict_target": artifact["verdict_target"],
        "verdict_branches": artifact["verdict_branches"],
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
