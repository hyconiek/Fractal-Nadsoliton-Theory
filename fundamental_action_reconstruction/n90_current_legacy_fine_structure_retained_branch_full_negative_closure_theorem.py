#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n90_current_legacy_fine_structure_retained_branch_full_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n86 = load_json(
        "fundamental_action_reconstruction/generated/n86_current_strict_side_literal_retention_obstruction_theorem_for_legacy_fine_structure_formula_summary.json"
    )
    n89 = load_json(
        "fundamental_action_reconstruction/generated/n89_current_strict_side_textual_retained_successor_obstruction_theorem_for_legacy_fine_structure_role_summary.json"
    )
    p85 = load_json(
        "fundamental_action_reconstruction/generated/p85_strict_side_object_lineage_upgrade_probe_for_legacy_fine_structure_role_summary.json"
    )

    checks_spec = [
        {
            "id": "literal_retention_closed",
            "actual": n86["theorem_result"]["literal_retention_path_closed_on_current_repo_state"],
            "expected": True,
            "meaning": "literal retention path is already closed negatively on the current repo state",
        },
        {
            "id": "textual_successor_closed",
            "actual": n89["theorem_result"]["textual_successor_path_closed_on_current_repo_state"],
            "expected": True,
            "meaning": "textual retained-successor path is already closed negatively on the current repo state",
        },
        {
            "id": "object_lineage_upgrade_absent",
            "actual": p85["object_lineage_upgrade_verdict_present"],
            "expected": False,
            "meaning": "object-lineage-upgrade verdict is absent on the current repo state",
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
            "step": "N90",
            "status": "N90_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_FINE_STRUCTURE_RETAINED_BRANCH_STATE",
            "scope": "current_legacy_fine_structure_retained_branch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N90",
            "status": "N90_DISCHARGED_CURRENT_LEGACY_FINE_STRUCTURE_RETAINED_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_fine_structure_retained_branch_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "retained_branch_closed_negatively_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branch": [
                "explicit_strict_side_replaced_verdict_for_the_legacy_fine_structure_role_by_an_explicit_strict_successor_semantics"
            ],
            "hard_limits": [
                "no_proof_that_the_replaced_branch_is_solved",
                "no_proof_that_the_retained_branch_is_impossible_forever",
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
