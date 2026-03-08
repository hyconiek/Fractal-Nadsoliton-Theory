#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "f57_first_post_additive_verdict_admissibility_branch_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n157 = load_json(
        "fundamental_action_reconstruction/generated/n157_current_additive_construction_attempt_verdict_layer_full_negative_closure_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n157_verdict_layer_closed",
            "actual": n157["theorem_result"][
                "additive_construction_attempt_verdict_layer_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N157 already closes the whole binary verdict layer negatively on the current repo state",
        },
        {
            "id": "n157_fixed_attempt_instance",
            "actual": n157["theorem_result"]["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N157 fixes the additive construction attempt under discussion",
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
            "stage": "F57",
            "lane": "post_additive_verdict_branch_order_only",
            "status": "F57_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_POST_ADDITIVE_VERDICT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F57",
            "lane": "post_additive_verdict_admissibility_first",
            "goal": "freeze_the_admissibility_branch_as_the_first_remaining_lower_branch_below_the_exhausted_additive_verdict_layer",
            "status": "F57_EXECUTED_FIRST_POST_ADDITIVE_VERDICT_ADMISSIBILITY_BRANCH_PACKET_NO_FALSE_PASS",
            "first_lower_branch_to_attack": "future_admissibility_test_of_a_future_constructed_source_object_for_S_sel_int_after_fixed_first_additive_attempt",
            "basis": [
                "binary_verdict_layer_already_closed_by_N157",
                "E_orient_still_presupposes_an_admissible_source_object",
                "downstream_B_sel_R_sel_O_sel_still_presupposes_source_and_export_stages",
                "preferred_order_remains_nadsoliton_to_light_to_matter_to_emergent_observer",
            ],
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT)


if __name__ == "__main__":
    main()
