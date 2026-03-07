#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f21_legacy_gravity_hierarchy_retained_branch_refinement_packet.json"
OUT_SUMMARY = GENERATED / "f21_legacy_gravity_hierarchy_retained_branch_refinement_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p93 = load_json(
        "fundamental_action_reconstruction/generated/p93_legacy_gravity_hierarchy_role_strict_side_branch_probe_summary.json"
    )

    literal_retention_present = False
    role_equivalence_retention_present = False

    checks_spec = [
        {
            "id": "p93_retained_branch_absent",
            "actual": p93["branch_state"]["retained_branch_present"],
            "expected": False,
            "meaning": "P93 already keeps the retained branch absent at branch level",
        },
        {
            "id": "literal_retention_present",
            "actual": literal_retention_present,
            "expected": False,
            "meaning": "the repo does not yet export literal retention of the old gravity-hierarchy claim on the strict side",
        },
        {
            "id": "role_equivalence_retention_present",
            "actual": role_equivalence_retention_present,
            "expected": False,
            "meaning": "the repo does not yet export role-equivalence retention for the gravity-hierarchy role on the strict side",
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
        "stage": "F21",
        "lane": "legacy_gravity_hierarchy_retained_branch_refinement_current_repo_state_only",
        "goal": "refine_the_missing_retained_branch_for_the_legacy_gravity_hierarchy_role_into_literal_vs_role_equivalence_retention",
        "status": "F21_EXECUTED_LEGACY_GRAVITY_HIERARCHY_RETAINED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "P93 already shows that the strict-side retained branch is absent, so the narrowest honest refinement is to split that missing retained verdict into literal retention vs role-equivalence retention, neither of which is currently exported",
        "retained_subbranch_state": {
            "literal_retention_present": literal_retention_present,
            "role_equivalence_retention_present": role_equivalence_retention_present,
        },
        "remaining_missing_objects": [
            "explicit_strict_side_literal_retention_of_exact_gravity_hierarchy_from_beta_to_the_N_scaling",
            "explicit_strict_side_role_equivalence_verdict_for_the_legacy_gravity_hierarchy_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F21",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "retained_subbranch_state": artifact["retained_subbranch_state"],
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
