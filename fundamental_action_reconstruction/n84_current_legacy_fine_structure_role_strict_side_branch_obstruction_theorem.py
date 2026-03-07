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
    / "n84_current_legacy_fine_structure_role_strict_side_branch_obstruction_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f13 = load_json(
        "fundamental_action_reconstruction/generated/f13_legacy_fine_structure_role_strict_side_partition_refinement_packet_summary.json"
    )
    p79 = load_json(
        "fundamental_action_reconstruction/generated/p79_legacy_fine_structure_role_strict_side_branch_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "f13_branch_refinement_present",
            "actual": f13["status"],
            "expected": "F13_EXECUTED_LEGACY_FINE_STRUCTURE_ROLE_STRICT_SIDE_PARTITION_REFINEMENT_PACKET_NO_FALSE_PASS",
            "meaning": "F13 refines the missing fine-structure verdict into retained and replaced branches",
        },
        {
            "id": "p79_probe_negative",
            "actual": p79["status"],
            "expected": "CURRENT_REPO_EXPORTS_NEITHER_RETAINED_NOR_REPLACED_STRICT_SIDE_BRANCH_FOR_THE_LEGACY_FINE_STRUCTURE_ROLE_AFTER_P79",
            "meaning": "P79 confirms that neither strict-side branch is currently exported",
        },
        {
            "id": "retained_branch_absent",
            "actual": p79["branch_state"]["retained_branch_present"],
            "expected": False,
            "meaning": "the strict-side retained branch is absent",
        },
        {
            "id": "replaced_branch_absent",
            "actual": p79["branch_state"]["replaced_branch_present"],
            "expected": False,
            "meaning": "the strict-side replaced branch is absent",
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
            "step": "N84",
            "status": "N84_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FINE_STRUCTURE_BRANCH_STATE",
            "scope": "current_legacy_fine_structure_role_branch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N84",
            "status": "N84_DISCHARGED_CURRENT_LEGACY_FINE_STRUCTURE_ROLE_STRICT_SIDE_BRANCH_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_fine_structure_role_branch_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "retained_branch_present": False,
                "replaced_branch_present": False,
                "strict_side_fine_structure_branch_verdict_is_unsupported": True,
                "full_closure_pass": False,
            },
            "missing_structure_classes": [
                "explicit_strict_side_retained_verdict_for_the_legacy_fine_structure_role",
                "explicit_strict_side_replaced_verdict_for_the_legacy_fine_structure_role_by_an_explicit_strict_successor_semantics",
            ],
            "hard_limits": [
                "no_proof_that_the_retained_branch_is_impossible_forever",
                "no_proof_that_the_replaced_branch_is_impossible_forever",
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
