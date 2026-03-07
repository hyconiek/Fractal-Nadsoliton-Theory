#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n112_current_strict_side_textual_method_successor_obstruction_theorem_for_legacy_gravity_hierarchy_role_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p105 = load_json(
        "fundamental_action_reconstruction/generated/p105_strict_side_textual_method_successor_probe_for_legacy_gravity_hierarchy_role_summary.json"
    )
    n110 = load_json(
        "fundamental_action_reconstruction/generated/n110_current_strict_side_object_successor_branch_full_negative_closure_theorem_for_legacy_gravity_hierarchy_role_summary.json"
    )

    checks_spec = [
        {
            "id": "n110_object_branch_closed",
            "actual": n110["theorem_result"]["object_successor_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "the full object-successor branch remains closed negatively on the current repo state",
        },
        {
            "id": "p105_textual_method_successor_verdict_absent",
            "actual": p105["textual_method_successor_verdict_present"],
            "expected": False,
            "meaning": "the textual method-successor sub-branch remains absent",
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
            "step": "N112",
            "status": "N112_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_SIDE_TEXTUAL_METHOD_SUCCESSOR_STATE",
            "scope": "current_strict_side_textual_method_successor_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N112",
            "status": "N112_DISCHARGED_CURRENT_STRICT_SIDE_TEXTUAL_METHOD_SUCCESSOR_OBSTRUCTION_THEOREM_FOR_LEGACY_GRAVITY_HIERARCHY_ROLE_NO_FALSE_PASS",
            "scope": "current_strict_side_textual_method_successor_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "object_successor_branch_closed_negatively_on_current_repo_state": True,
                "textual_method_successor_subbranch_present": False,
                "method_successor_frontier_still_open": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "explicit_method_lineage_upgrade_verdict_elevating_the_existing_qw2115_micro_supported_beta_hierarchy_bridge_chain_into_replacement_semantics_for_the_legacy_gravity_hierarchy_role"
            ],
            "hard_limits": [
                "no_proof_that_the_method_lineage_upgrade_subbranch_is_impossible_forever",
                "no_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT)


if __name__ == "__main__":
    main()
