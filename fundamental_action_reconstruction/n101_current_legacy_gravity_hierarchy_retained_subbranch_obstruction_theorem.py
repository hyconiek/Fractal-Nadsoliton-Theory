#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n101_current_legacy_gravity_hierarchy_retained_subbranch_obstruction_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f21 = load_json(
        "fundamental_action_reconstruction/generated/f21_legacy_gravity_hierarchy_retained_branch_refinement_packet_summary.json"
    )
    p94 = load_json(
        "fundamental_action_reconstruction/generated/p94_legacy_gravity_hierarchy_retained_subbranch_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "f21_retained_subbranch_refinement_present",
            "actual": f21["status"],
            "expected": "F21_EXECUTED_LEGACY_GRAVITY_HIERARCHY_RETAINED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
            "meaning": "F21 refines the retained branch into literal-retention and role-equivalence subbranches",
        },
        {
            "id": "p94_probe_negative",
            "actual": p94["status"],
            "expected": "CURRENT_REPO_EXPORTS_NEITHER_LITERAL_NOR_ROLE_EQUIVALENCE_RETAINED_SUBBRANCH_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P94",
            "meaning": "P94 confirms that neither retained sub-branch is currently exported",
        },
        {
            "id": "literal_retention_absent",
            "actual": p94["retained_subbranch_state"]["literal_retention_present"],
            "expected": False,
            "meaning": "literal retention is absent",
        },
        {
            "id": "role_equivalence_retention_absent",
            "actual": p94["retained_subbranch_state"]["role_equivalence_retention_present"],
            "expected": False,
            "meaning": "role-equivalence retention is absent",
        },
    ]

    checks: list[dict[str, Any]] = []
    mismatches: list[str] = []
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
            "step": "N101",
            "status": "N101_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_GRAVITY_HIERARCHY_RETAINED_SUBBRANCH_STATE",
            "scope": "current_legacy_gravity_hierarchy_retained_subbranch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N101",
            "status": "N101_DISCHARGED_CURRENT_LEGACY_GRAVITY_HIERARCHY_RETAINED_SUBBRANCH_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_gravity_hierarchy_retained_subbranch_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "literal_retention_present": False,
                "role_equivalence_retention_present": False,
                "retained_subbranch_is_unsupported": True,
                "full_closure_pass": False,
            },
            "missing_structure_classes": [
                "explicit_strict_side_literal_retention_of_exact_gravity_hierarchy_from_beta_to_the_N_scaling",
                "explicit_strict_side_role_equivalence_verdict_for_the_legacy_gravity_hierarchy_role",
            ],
            "hard_limits": [
                "no_proof_that_literal_retention_is_impossible_forever",
                "no_proof_that_role_equivalence_retention_is_impossible_forever",
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
