#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n77_current_strict_side_textual_object_successor_obstruction_theorem_for_legacy_weinberg_role_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p74 = load_json(
        "fundamental_action_reconstruction/generated/p74_strict_side_textual_object_successor_probe_for_legacy_weinberg_role_summary.json"
    )
    n73 = load_json(
        "fundamental_action_reconstruction/generated/n73_current_legacy_weinberg_retained_branch_full_negative_closure_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n73_retained_branch_closed",
            "actual": n73["theorem_result"]["retained_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "the retained branch remains closed negatively on the current repo state",
        },
        {
            "id": "p74_textual_object_successor_verdict_absent",
            "actual": p74["textual_object_successor_verdict_present"],
            "expected": False,
            "meaning": "the textual object-successor sub-branch remains absent",
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
            "step": "N77",
            "status": "N77_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_SIDE_TEXTUAL_OBJECT_SUCCESSOR_STATE",
            "scope": "current_strict_side_textual_object_successor_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N77",
            "status": "N77_DISCHARGED_CURRENT_STRICT_SIDE_TEXTUAL_OBJECT_SUCCESSOR_OBSTRUCTION_THEOREM_FOR_LEGACY_WEINBERG_ROLE_NO_FALSE_PASS",
            "scope": "current_strict_side_textual_object_successor_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "retained_branch_closed_negatively_on_current_repo_state": True,
                "textual_object_successor_subbranch_present": False,
                "object_successor_frontier_still_open": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "explicit_object_lineage_upgrade_verdict_elevating_the_existing_sin2_theta_w_mz_candidate_chain_into_replacement_semantics_for_the_legacy_weinberg_angle_role"
            ],
            "method_branch_still_open": [
                "explicit_method_successor_semantics_verdict_identifying_qw2098_sin2_from_nonanchor_ew_pole_chain_as_the_strict_side_successor_semantics_replacing_the_legacy_weinberg_angle_role"
            ],
            "hard_limits": [
                "no_proof_that_the_object_lineage_upgrade_subbranch_is_impossible_forever",
                "no_method_successor_semantics_closure",
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
