#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "f58_first_post_additive_admissibility_orientation_export_branch_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n158 = load_json(
        "fundamental_action_reconstruction/generated/n158_current_additive_post_verdict_admissibility_branch_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n158_additive_admissibility_branch_closed",
            "actual": n158["theorem_result"]["current_additive_post_verdict_admissibility_branch_obstructed"],
            "expected": True,
            "meaning": "N158 already closes the additive-specific admissibility branch negatively on the current repo state",
        },
        {
            "id": "n158_fixed_attempt_instance",
            "actual": n158["theorem_result"]["fixed_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "N158 fixes the additive construction attempt under discussion",
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
            "stage": "F58",
            "lane": "post_additive_admissibility_branch_order_only",
            "status": "F58_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_POST_ADDITIVE_ADMISSIBILITY_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F58",
            "lane": "post_additive_admissibility_orientation_first",
            "goal": "freeze_the_additive_specific_E_orient_branch_as_the_first_remaining_lower_branch_after_additive_admissibility_obstruction",
            "status": "F58_EXECUTED_FIRST_POST_ADDITIVE_ADMISSIBILITY_ORIENTATION_EXPORT_BRANCH_PACKET_NO_FALSE_PASS",
            "first_remaining_lower_branch_to_attack": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object_for_the_fixed_first_additive_attempt",
            "basis": [
                "additive_specific_admissibility_branch_already_closed_by_N158",
                "F32_already_freezes_the_admissible_contract_for_E_orient",
                "downstream_B_sel_R_sel_O_sel_still_presupposes_source_and_orientation_export",
                "preferred_order_remains_nadsoliton_to_light_to_matter_to_emergent_observer",
                "observer_information_deficit_cannot_be_promoted_above_E_orient",
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
