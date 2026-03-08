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
    / "f54_first_future_additive_construction_attempt_verdict_branch_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f54_first_future_additive_construction_attempt_verdict_branch_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f53 = load_json(
        "fundamental_action_reconstruction/generated/f53_first_future_additive_construction_attempt_verdict_target_packet_summary.json"
    )
    n153 = load_json(
        "fundamental_action_reconstruction/generated/n153_next_constructive_move_reduced_to_one_first_additive_construction_attempt_verdict_target_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n153_verdict_target",
            "actual": n153["theorem_result"]["explicit_future_additive_attempt_verdict_target"],
            "expected": "success_or_failure_verdict(construct_attempt_v1(S_sel_int_additive_attempt_target_v1))",
            "meaning": "N153 already fixes one additive-construction verdict target",
        },
        {
            "id": "f53_success_verdict_absent",
            "actual": f53["future_additive_attempt_verdict_target"]["success_verdict_present"],
            "expected": False,
            "meaning": "F53 does not overclaim a success verdict",
        },
        {
            "id": "f53_failure_verdict_absent",
            "actual": f53["future_additive_attempt_verdict_target"]["failure_verdict_present"],
            "expected": False,
            "meaning": "F53 does not overclaim a failure verdict",
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
        "stage": "F54",
        "lane": "future_additive_attempt_verdict_branch_refinement_only",
        "goal": "freeze_the_two_explicit_verdict_branches_on_the_fixed_additive_construction_attempt_verdict_target",
        "status": "F54_EXECUTED_FIRST_FUTURE_ADDITIVE_CONSTRUCTION_ATTEMPT_VERDICT_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
        "verdict_branch_split": {
            "fixed_verdict_target": "success_or_failure_verdict(construct_attempt_v1(S_sel_int_additive_attempt_target_v1))",
            "success_branch": "explicit_success_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "failure_branch": "explicit_failure_verdict_for_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F54",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "verdict_branch_split": artifact["verdict_branch_split"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
