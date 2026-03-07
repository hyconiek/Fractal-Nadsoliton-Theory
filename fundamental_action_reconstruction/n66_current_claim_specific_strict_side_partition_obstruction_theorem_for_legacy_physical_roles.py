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
    / "n66_current_claim_specific_strict_side_partition_obstruction_theorem_for_legacy_physical_roles_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f5 = load_json(
        "fundamental_action_reconstruction/generated/f5_legacy_physical_role_partition_refinement_packet_summary.json"
    )
    p63 = load_json(
        "fundamental_action_reconstruction/generated/p63_claim_specific_strict_side_partition_probe_for_legacy_physical_roles_summary.json"
    )

    checks_spec = [
        {
            "id": "f5_packet_present",
            "actual": f5["status"],
            "expected": "F5_EXECUTED_LEGACY_PHYSICAL_ROLE_PARTITION_REFINEMENT_PACKET_NO_FALSE_PASS",
            "meaning": "F5 refines the broad partition blocker into three claim-specific verdict gaps",
        },
        {
            "id": "p63_probe_negative",
            "actual": p63["status"],
            "expected": "CURRENT_REPO_EXPORTS_BROAD_RETAINED_VS_REPLACED_PARTITION_BUT_NO_CLAIM_SPECIFIC_STRICT_SIDE_PARTITION_FOR_LEGACY_PHYSICAL_ROLES_AFTER_P63",
            "meaning": "P63 confirms that only the broad partition exists while claim-specific strict-side verdicts remain absent",
        },
        {
            "id": "weinberg_partition_verdict_absent",
            "actual": p63["claim_specific_partition_state"][
                "legacy_weinberg_angle_role_partition_verdict_present"
            ],
            "expected": False,
            "meaning": "no strict-side claim-specific partition verdict exists for the legacy Weinberg-angle role",
        },
        {
            "id": "fine_structure_partition_verdict_absent",
            "actual": p63["claim_specific_partition_state"][
                "legacy_fine_structure_role_partition_verdict_present"
            ],
            "expected": False,
            "meaning": "no strict-side claim-specific partition verdict exists for the legacy fine-structure role",
        },
        {
            "id": "gravity_hierarchy_partition_verdict_absent",
            "actual": p63["claim_specific_partition_state"][
                "legacy_gravity_hierarchy_role_partition_verdict_present"
            ],
            "expected": False,
            "meaning": "no strict-side claim-specific partition verdict exists for the legacy gravity-hierarchy role",
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
            "step": "N66",
            "status": "N66_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CLAIM_SPECIFIC_STRICT_SIDE_PARTITION_STATE",
            "scope": "current_claim_specific_strict_side_partition_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N66",
            "status": "N66_DISCHARGED_CURRENT_CLAIM_SPECIFIC_STRICT_SIDE_PARTITION_OBSTRUCTION_THEOREM_FOR_LEGACY_PHYSICAL_ROLES_NO_FALSE_PASS",
            "scope": "current_claim_specific_strict_side_partition_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "only_broad_partition_present": True,
                "claim_specific_weinberg_partition_verdict_present": False,
                "claim_specific_fine_structure_partition_verdict_present": False,
                "claim_specific_gravity_hierarchy_partition_verdict_present": False,
                "claim_specific_strict_side_partition_is_unsupported": True,
                "full_closure_pass": False,
            },
            "missing_structure_classes": [
                "explicit_strict_side_retained_or_replaced_verdict_for_the_legacy_weinberg_angle_role",
                "explicit_strict_side_retained_or_replaced_verdict_for_the_legacy_fine_structure_role",
                "explicit_strict_side_retained_or_replaced_verdict_for_the_legacy_gravity_hierarchy_role",
            ],
            "hard_limits": [
                "no_proof_that_any_role_is_retained",
                "no_proof_that_any_role_is_replaced",
                "no_proof_that_no_future_claim_specific_partition_can_ever_exist",
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
