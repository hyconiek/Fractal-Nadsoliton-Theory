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
    / "f20_legacy_gravity_hierarchy_role_strict_side_partition_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f20_legacy_gravity_hierarchy_role_strict_side_partition_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p62 = load_json(
        "fundamental_action_reconstruction/generated/p62_legacy_physical_role_transfer_to_strict_gate_kernel_probe_summary.json"
    )
    p63 = load_json(
        "fundamental_action_reconstruction/generated/p63_claim_specific_strict_side_partition_probe_for_legacy_physical_roles_summary.json"
    )

    retained_branch_present = False
    replaced_branch_present = False

    checks_spec = [
        {
            "id": "p62_gravity_hierarchy_role_transfer_absent",
            "actual": p62["role_transfer_state"]["gravity_hierarchy_role_transfer_present"],
            "expected": False,
            "meaning": "the retained branch is not silently discharged through role transfer onto K_strict_gate",
        },
        {
            "id": "p63_gravity_hierarchy_claim_specific_partition_verdict_absent",
            "actual": p63["claim_specific_partition_state"][
                "legacy_gravity_hierarchy_role_partition_verdict_present"
            ],
            "expected": False,
            "meaning": "the repo still exports no claim-specific strict-side verdict for the gravity-hierarchy role",
        },
        {
            "id": "retained_branch_present",
            "actual": retained_branch_present,
            "expected": False,
            "meaning": "no explicit strict-side retained verdict for the legacy gravity-hierarchy role is currently exported",
        },
        {
            "id": "replaced_branch_present",
            "actual": replaced_branch_present,
            "expected": False,
            "meaning": "no explicit strict-side replaced verdict for the legacy gravity-hierarchy role is currently exported",
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
        "stage": "F20",
        "lane": "legacy_gravity_hierarchy_role_strict_side_partition_refinement_current_repo_state_only",
        "goal": "refine_the_missing_strict_side_gravity_hierarchy_verdict_into_retained_vs_replaced_branches",
        "status": "F20_EXECUTED_LEGACY_GRAVITY_HIERARCHY_ROLE_STRICT_SIDE_PARTITION_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "P62 keeps the gravity-hierarchy role-transfer branch absent and P63 keeps the claim-specific strict-side verdict absent, so the narrowest honest refinement is to split that one missing verdict into a retained branch and a replaced branch, neither of which is currently exported",
        "branch_state": {
            "retained_branch_present": retained_branch_present,
            "replaced_branch_present": replaced_branch_present,
        },
        "remaining_missing_objects": [
            "explicit_strict_side_retained_verdict_for_the_legacy_gravity_hierarchy_role",
            "explicit_strict_side_replaced_verdict_for_the_legacy_gravity_hierarchy_role_by_an_explicit_strict_successor_semantics",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F20",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "branch_state": artifact["branch_state"],
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
