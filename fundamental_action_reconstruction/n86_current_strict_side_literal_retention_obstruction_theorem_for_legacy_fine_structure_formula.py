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
    / "n86_current_strict_side_literal_retention_obstruction_theorem_for_legacy_fine_structure_formula_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p81 = load_json(
        "fundamental_action_reconstruction/generated/p81_strict_side_literal_retention_probe_for_legacy_fine_structure_formula_summary.json"
    )

    checks_spec = [
        {
            "id": "p81_probe_negative",
            "actual": p81["status"],
            "expected": "CURRENT_STRICT_SIDE_AUTHORITATIVE_SOURCES_DO_NOT_EXPORT_LITERAL_RETENTION_OF_THE_LEGACY_FINE_STRUCTURE_FORMULA_AFTER_P81",
            "meaning": "P81 confirms that strict-side authoritative sources do not literally retain the old fine-structure formula",
        },
        {
            "id": "literal_retention_absent",
            "actual": p81["literal_retention_present"],
            "expected": False,
            "meaning": "literal strict-side retention is absent on the current repo state",
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
            "step": "N86",
            "status": "N86_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LITERAL_RETENTION_STATE",
            "scope": "current_strict_side_literal_retention_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N86",
            "status": "N86_DISCHARGED_CURRENT_STRICT_SIDE_LITERAL_RETENTION_OBSTRUCTION_THEOREM_FOR_LEGACY_FINE_STRUCTURE_FORMULA_NO_FALSE_PASS",
            "scope": "current_strict_side_literal_retention_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "literal_retention_present": False,
                "literal_retention_path_closed_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branch": [
                "explicit_strict_side_role_equivalence_verdict_for_the_legacy_fine_structure_role"
            ],
            "hard_limits": [
                "no_proof_that_role_equivalence_retention_is_absent",
                "no_proof_that_the_replaced_branch_is_solved",
                "no_proof_that_literal_retention_is_impossible_forever",
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
