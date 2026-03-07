#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p79_legacy_fine_structure_role_strict_side_branch_probe.json"
OUT_SUMMARY = (
    GENERATED / "p79_legacy_fine_structure_role_strict_side_branch_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f13 = load_json(
        "fundamental_action_reconstruction/generated/f13_legacy_fine_structure_role_strict_side_partition_refinement_packet_summary.json"
    )

    branch_state = f13["branch_state"]

    checks_spec = [
        {
            "id": "f13_partition_refinement_present",
            "actual": f13["status"],
            "expected": "F13_EXECUTED_LEGACY_FINE_STRUCTURE_ROLE_STRICT_SIDE_PARTITION_REFINEMENT_PACKET_NO_FALSE_PASS",
            "meaning": "F13 refines the missing fine-structure verdict into retained vs replaced branches",
        },
        {
            "id": "retained_branch_present",
            "actual": branch_state["retained_branch_present"],
            "expected": False,
            "meaning": "the repo does not yet export the retained branch for the legacy fine-structure role",
        },
        {
            "id": "replaced_branch_present",
            "actual": branch_state["replaced_branch_present"],
            "expected": False,
            "meaning": "the repo does not yet export the replaced branch for the legacy fine-structure role",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "P79",
        "lane": "legacy_fine_structure_role_strict_side_branch_probe_current_repo_state_only",
        "goal": "test_whether_the_current_repo_already_exports_either_the_retained_or_the_replaced_strict_side_branch_for_the_legacy_fine_structure_role",
        "status": "CURRENT_REPO_EXPORTS_NEITHER_RETAINED_NOR_REPLACED_STRICT_SIDE_BRANCH_FOR_THE_LEGACY_FINE_STRUCTURE_ROLE_AFTER_P79",
        "reason": "F13 already reduces the missing claim-specific fine-structure verdict to two branches, and the current repo exports neither one",
        "branch_state": branch_state,
        "remaining_missing_objects": [
            "explicit_strict_side_retained_verdict_for_the_legacy_fine_structure_role",
            "explicit_strict_side_replaced_verdict_for_the_legacy_fine_structure_role_by_an_explicit_strict_successor_semantics",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P79",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "branch_state": branch_state,
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
