#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p65_legacy_weinberg_retained_subbranch_probe.json"
OUT_SUMMARY = GENERATED / "p65_legacy_weinberg_retained_subbranch_probe_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f7 = load_json(
        "fundamental_action_reconstruction/generated/f7_legacy_weinberg_retained_branch_refinement_packet_summary.json"
    )

    subbranch_state = f7["retained_subbranch_state"]

    checks_spec = [
        {
            "id": "f7_retained_refinement_present",
            "actual": f7["status"],
            "expected": "F7_EXECUTED_LEGACY_WEINBERG_RETAINED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
            "meaning": "F7 refines the retained branch into literal vs role-equivalence retention",
        },
        {
            "id": "literal_retention_present",
            "actual": subbranch_state["literal_retention_present"],
            "expected": False,
            "meaning": "the repo does not yet export literal retention of the old Weinberg formula on the strict side",
        },
        {
            "id": "role_equivalence_retention_present",
            "actual": subbranch_state["role_equivalence_retention_present"],
            "expected": False,
            "meaning": "the repo does not yet export role-equivalence retention for the Weinberg role on the strict side",
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
        "stage": "P65",
        "lane": "legacy_weinberg_retained_subbranch_probe_current_repo_state_only",
        "goal": "test_whether_the_current_repo_already_exports_either_literal_or_role_equivalence_retention_for_the_legacy_weinberg_angle_role",
        "status": "CURRENT_REPO_EXPORTS_NEITHER_LITERAL_NOR_ROLE_EQUIVALENCE_RETAINED_SUBBRANCH_FOR_THE_LEGACY_WEINBERG_ANGLE_ROLE_AFTER_P65",
        "reason": "F7 already reduces the retained branch to literal-retention and role-equivalence-retention subbranches, and the current repo exports neither one",
        "retained_subbranch_state": subbranch_state,
        "remaining_missing_objects": [
            "explicit_strict_side_literal_retention_of_sin2_thetaW_equals_alpha_geo_over_12",
            "explicit_strict_side_role_equivalence_verdict_for_the_legacy_weinberg_angle_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P65",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "retained_subbranch_state": subbranch_state,
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
