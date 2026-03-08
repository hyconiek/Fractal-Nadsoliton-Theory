#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f71_first_post_contract_compliant_admissibility_orientation_export_branch_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f32 = load_json(
        "fundamental_action_reconstruction/generated/f32_initial_future_strict_core_orientation_export_admission_packet_summary.json"
    )
    n173 = load_json(
        "fundamental_action_reconstruction/generated/n173_current_contract_compliant_post_verdict_admissibility_branch_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n173_admissibility_obstructed",
            "actual": n173["theorem_result"]["current_contract_compliant_post_verdict_admissibility_branch_obstructed"],
            "expected": True,
            "meaning": "N173 already obstructs the contract-compliant admissibility branch",
        },
        {
            "id": "f32_orientation_export_name",
            "actual": f32["orientation_export_name"],
            "expected": "E_orient",
            "meaning": "F32 already fixes the orientation export node name",
        },
        {
            "id": "f32_upstream_source_name",
            "actual": f32["upstream_source_name"],
            "expected": "S_sel_int",
            "meaning": "F32 already fixes the upstream source node name",
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
            "stage": "F71",
            "lane": "post_contract_compliant_admissibility_branch_order_only",
            "status": "F71_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_POST_CONTRACT_COMPLIANT_ADMISSIBILITY_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F71",
            "lane": "post_contract_compliant_admissibility_orientation_export_first",
            "goal": "freeze_the_orientation_export_branch_as_the_first_remaining_lower_branch_after_the_contract_compliant_admissibility_obstruction",
            "status": "F71_EXECUTED_FIRST_POST_CONTRACT_COMPLIANT_ADMISSIBILITY_ORIENTATION_EXPORT_BRANCH_PACKET_NO_FALSE_PASS",
            "fixed_attempt_instance": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "first_remaining_lower_branch": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object_for_the_fixed_first_contract_compliant_additive_attempt",
            "basis": [
                "contract_compliant_admissibility_branch_already_obstructed_by_N173",
                "F32_already_freezes_the_E_orient_admission_contract",
                "downstream_B_sel_R_sel_O_sel_still_presupposes_source_and_orientation_export",
                "preferred_order_remains_nadsoliton_to_light_to_matter_to_emergent_observer",
            ],
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
