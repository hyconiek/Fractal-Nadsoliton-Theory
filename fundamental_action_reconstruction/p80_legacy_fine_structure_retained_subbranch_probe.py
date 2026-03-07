#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p80_legacy_fine_structure_retained_subbranch_probe.json"
OUT_SUMMARY = GENERATED / "p80_legacy_fine_structure_retained_subbranch_probe_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f14 = load_json(
        "fundamental_action_reconstruction/generated/f14_legacy_fine_structure_retained_branch_refinement_packet_summary.json"
    )

    subbranch_state = f14["retained_subbranch_state"]

    checks_spec = [
        {
            "id": "f14_retained_refinement_present",
            "actual": f14["status"],
            "expected": "F14_EXECUTED_LEGACY_FINE_STRUCTURE_RETAINED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
            "meaning": "F14 refines the retained branch into literal vs role-equivalence retention",
        },
        {
            "id": "literal_retention_present",
            "actual": subbranch_state["literal_retention_present"],
            "expected": False,
            "meaning": "the repo does not yet export literal retention of the old fine-structure formula on the strict side",
        },
        {
            "id": "role_equivalence_retention_present",
            "actual": subbranch_state["role_equivalence_retention_present"],
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
        "stage": "P80",
        "lane": "legacy_fine_structure_retained_subbranch_probe_current_repo_state_only",
        "goal": "test_whether_the_current_repo_already_exports_either_literal_or_role_equivalence_retention_for_the_legacy_fine_structure_role",
        "status": "CURRENT_REPO_EXPORTS_NEITHER_LITERAL_NOR_ROLE_EQUIVALENCE_RETAINED_SUBBRANCH_FOR_THE_LEGACY_FINE_STRUCTURE_ROLE_AFTER_P80",
        "reason": "F14 already reduces the retained branch to literal-retention and role-equivalence-retention subbranches, and the current repo exports neither one",
        "retained_subbranch_state": subbranch_state,
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
        "stage": "P80",
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
