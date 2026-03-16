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
    / "f643_first_post_verdict_admissibility_branch_packet_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f643_first_post_verdict_admissibility_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n532 = load_json(
        "fundamental_action_reconstruction/generated/n532_current_failure_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
    )
    n533 = load_json(
        "fundamental_action_reconstruction/generated/n533_current_success_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "n532_failure_branch_obstruction_discharged",
            "actual": n532["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N532 already packages the current failure-side obstruction (v1)",
        },
        {
            "id": "n533_success_branch_obstruction_discharged",
            "actual": n533["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N533 already packages the current success-side obstruction (v1)",
        },
        {
            "id": "n533_first_remaining_open_branch",
            "actual": n533["remaining_open_branches"][0],
            "expected": "future_admissibility_test_of_a_future_constructed_source_object",
            "meaning": "after the verdict layer, the first remaining lower branch is admissibility (v1)",
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
        "stage": "F643",
        "lane": "post_verdict_admissibility_branch_first_only",
        "goal": "freeze_the_admissibility_branch_as_the_first_remaining_lower_branch_below_the_exhausted_binary_verdict_layer_for_seed_v1",
        "status": "F643_EXECUTED_FIRST_POST_VERDICT_ADMISSIBILITY_BRANCH_PACKET_FOR_SEED_V1_NO_FALSE_PASS",
        "first_lower_branch_to_attack": "future_admissibility_test_of_a_future_constructed_source_object",
        "branch_ordering_basis": [
            "failure_side_obstruction_already_packaged_by_N532",
            "success_side_obstruction_already_packaged_by_N533",
            "E_orient_requires_an_admissible_source_object",
            "downstream_B_sel_R_sel_O_sel_requires_source_object_and_export_stage",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F643",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "first_lower_branch_to_attack": artifact["first_lower_branch_to_attack"],
        "branch_ordering_basis": artifact["branch_ordering_basis"],
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

