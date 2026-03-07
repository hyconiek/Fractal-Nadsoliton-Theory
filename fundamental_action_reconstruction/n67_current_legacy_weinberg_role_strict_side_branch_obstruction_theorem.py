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
    / "n67_current_legacy_weinberg_role_strict_side_branch_obstruction_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f6 = load_json(
        "fundamental_action_reconstruction/generated/f6_legacy_weinberg_role_strict_side_partition_refinement_packet_summary.json"
    )
    p64 = load_json(
        "fundamental_action_reconstruction/generated/p64_legacy_weinberg_role_strict_side_branch_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "f6_branch_refinement_present",
            "actual": f6["status"],
            "expected": "F6_EXECUTED_LEGACY_WEINBERG_ROLE_STRICT_SIDE_PARTITION_REFINEMENT_PACKET_NO_FALSE_PASS",
            "meaning": "F6 refines the missing Weinberg verdict into retained and replaced branches",
        },
        {
            "id": "p64_probe_negative",
            "actual": p64["status"],
            "expected": "CURRENT_REPO_EXPORTS_NEITHER_RETAINED_NOR_REPLACED_STRICT_SIDE_BRANCH_FOR_THE_LEGACY_WEINBERG_ANGLE_ROLE_AFTER_P64",
            "meaning": "P64 confirms that neither strict-side branch is currently exported",
        },
        {
            "id": "retained_branch_absent",
            "actual": p64["branch_state"]["retained_branch_present"],
            "expected": False,
            "meaning": "the strict-side retained branch is absent",
        },
        {
            "id": "replaced_branch_absent",
            "actual": p64["branch_state"]["replaced_branch_present"],
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
            "step": "N67",
            "status": "N67_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_WEINBERG_BRANCH_STATE",
            "scope": "current_legacy_weinberg_role_branch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N67",
            "status": "N67_DISCHARGED_CURRENT_LEGACY_WEINBERG_ROLE_STRICT_SIDE_BRANCH_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_weinberg_role_branch_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "retained_branch_present": False,
                "replaced_branch_present": False,
                "strict_side_weinberg_branch_verdict_is_unsupported": True,
                "full_closure_pass": False,
            },
            "missing_structure_classes": [
                "explicit_strict_side_retained_verdict_for_the_legacy_weinberg_angle_role",
                "explicit_strict_side_replaced_verdict_for_the_legacy_weinberg_angle_role_by_an_explicit_strict_successor_semantics",
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
