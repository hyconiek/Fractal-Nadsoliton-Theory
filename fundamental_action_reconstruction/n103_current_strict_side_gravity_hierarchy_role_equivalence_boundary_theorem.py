#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n103_current_strict_side_gravity_hierarchy_role_equivalence_boundary_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p96 = load_json(
        "fundamental_action_reconstruction/generated/p96_strict_side_role_equivalence_probe_for_legacy_gravity_hierarchy_role_summary.json"
    )

    checks_spec = [
        {
            "id": "p96_status_matches_expected_boundary",
            "actual": p96["status"],
            "expected": "CURRENT_REPO_EXPORTS_STRICT_SIDE_GRAVITY_HIERARCHY_BETA20_CANDIDATE_BUT_NO_EXPLICIT_LEGACY_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_VERDICT_AFTER_P96",
            "meaning": "P96 confirms that the strict side exports a candidate object but no explicit role-equivalence verdict",
        },
        {
            "id": "strict_side_candidate_object_present",
            "actual": p96["strict_side_candidate_object"]["present"],
            "expected": True,
            "meaning": "the strict-side gravity-hierarchy candidate object is really exported",
        },
        {
            "id": "explicit_role_equivalence_verdict_absent",
            "actual": p96["explicit_role_equivalence_verdict_present"],
            "expected": False,
            "meaning": "the explicit legacy-to-strict role-equivalence verdict is still absent",
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
            "step": "N103",
            "status": "N103_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_SIDE_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_STATE",
            "scope": "current_strict_side_gravity_hierarchy_role_equivalence_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N103",
            "status": "N103_DISCHARGED_CURRENT_STRICT_SIDE_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_BOUNDARY_THEOREM_NO_FALSE_PASS",
            "scope": "current_strict_side_gravity_hierarchy_role_equivalence_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "strict_side_candidate_object_present": True,
                "explicit_role_equivalence_verdict_present": False,
                "retained_branch_fully_discharged": False,
                "full_closure_pass": False,
            },
            "remaining_open_branch": [
                "explicit_legacy_to_strict_semantic_transfer_verdict_identifying_gravity_hierarchy_beta20_as_the_retained_strict_side_successor_of_the_legacy_gravity_hierarchy_role"
            ],
            "hard_limits": [
                "no_proof_that_the_retained_branch_is_impossible_forever",
                "no_proof_that_the_replaced_branch_is_solved",
                "no_proof_that_gravity_hierarchy_beta20_cannot_become_the_retained_successor_later",
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
