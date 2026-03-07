#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f14_legacy_fine_structure_retained_branch_refinement_packet.json"
OUT_SUMMARY = GENERATED / "f14_legacy_fine_structure_retained_branch_refinement_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p79 = load_json(
        "fundamental_action_reconstruction/generated/p79_legacy_fine_structure_role_strict_side_branch_probe_summary.json"
    )

    literal_retention_present = False
    role_equivalence_retention_present = False

    checks_spec = [
        {
            "id": "p79_retained_branch_absent",
            "actual": p79["branch_state"]["retained_branch_present"],
            "expected": False,
            "meaning": "P79 already keeps the retained branch absent at branch level",
        },
        {
            "id": "literal_retention_present",
            "actual": literal_retention_present,
            "expected": False,
            "meaning": "the repo does not yet export literal retention of the old fine-structure formula on the strict side",
        },
        {
            "id": "role_equivalence_retention_present",
            "actual": role_equivalence_retention_present,
            "expected": False,
            "meaning": "the repo does not yet export role-equivalence retention for the fine-structure role on the strict side",
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
        "stage": "F14",
        "lane": "legacy_fine_structure_retained_branch_refinement_current_repo_state_only",
        "goal": "refine_the_missing_retained_branch_for_the_legacy_fine_structure_role_into_literal_vs_role_equivalence_retention",
        "status": "F14_EXECUTED_LEGACY_FINE_STRUCTURE_RETAINED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "P79 already shows that the strict-side retained branch is absent, so the narrowest honest refinement is to split that missing retained verdict into literal retention vs role-equivalence retention, neither of which is currently exported",
        "retained_subbranch_state": {
            "literal_retention_present": literal_retention_present,
            "role_equivalence_retention_present": role_equivalence_retention_present,
        },
        "remaining_missing_objects": [
            "explicit_strict_side_literal_retention_of_alpha_em_inverse_equals_alpha_geo_over_2beta_tors_times_1_minus_beta_tors",
            "explicit_strict_side_role_equivalence_verdict_for_the_legacy_fine_structure_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F14",
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
