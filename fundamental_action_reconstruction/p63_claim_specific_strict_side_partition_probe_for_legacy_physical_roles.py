#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p63_claim_specific_strict_side_partition_probe_for_legacy_physical_roles.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p63_claim_specific_strict_side_partition_probe_for_legacy_physical_roles_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f5 = load_json(
        "fundamental_action_reconstruction/generated/f5_legacy_physical_role_partition_refinement_packet_summary.json"
    )

    partition_state = f5["claim_specific_partition_state"]

    checks_spec = [
        {
            "id": "f5_partition_refinement_packet_present",
            "actual": f5["status"],
            "expected": "F5_EXECUTED_LEGACY_PHYSICAL_ROLE_PARTITION_REFINEMENT_PACKET_NO_FALSE_PASS",
            "meaning": "F5 reduces the broad partition blocker to three claim-specific verdict gaps",
        },
        {
            "id": "weinberg_partition_verdict_present",
            "actual": partition_state[
                "legacy_weinberg_angle_role_partition_verdict_present"
            ],
            "expected": False,
            "meaning": "the repo does not yet export a strict-side partition verdict for the legacy Weinberg-angle role",
        },
        {
            "id": "fine_structure_partition_verdict_present",
            "actual": partition_state[
                "legacy_fine_structure_role_partition_verdict_present"
            ],
            "expected": False,
            "meaning": "the repo does not yet export a strict-side partition verdict for the legacy fine-structure role",
        },
        {
            "id": "gravity_hierarchy_partition_verdict_present",
            "actual": partition_state[
                "legacy_gravity_hierarchy_role_partition_verdict_present"
            ],
            "expected": False,
            "meaning": "the repo does not yet export a strict-side partition verdict for the legacy gravity-hierarchy role",
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

    remaining_missing = [
        "explicit_strict_side_retained_or_replaced_verdict_for_the_legacy_weinberg_angle_role",
        "explicit_strict_side_retained_or_replaced_verdict_for_the_legacy_fine_structure_role",
        "explicit_strict_side_retained_or_replaced_verdict_for_the_legacy_gravity_hierarchy_role",
    ]

    artifact = {
        "stage": "P63",
        "lane": "claim_specific_strict_side_partition_probe_for_legacy_physical_roles_current_repo_state_only",
        "goal": "test_whether_the_current_repo_already_exports_claim_specific_strict_side_retained_or_replaced_verdicts_for_the_three_legacy_physical_roles",
        "status": "CURRENT_REPO_EXPORTS_BROAD_RETAINED_VS_REPLACED_PARTITION_BUT_NO_CLAIM_SPECIFIC_STRICT_SIDE_PARTITION_FOR_LEGACY_PHYSICAL_ROLES_AFTER_P63",
        "reason": "F5 already reduces the broad partition blocker to three claim-specific verdict gaps, and the current repo still exports none of those three strict-side verdicts",
        "claim_specific_partition_state": partition_state,
        "remaining_missing_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P63",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "claim_specific_partition_state": partition_state,
        "remaining_missing_objects": remaining_missing,
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
