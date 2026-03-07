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
    / "n65_current_legacy_physical_role_transfer_to_strict_gate_kernel_obstruction_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f4 = load_json(
        "fundamental_action_reconstruction/generated/f4_legacy_physical_role_transfer_classification_packet_summary.json"
    )
    p62 = load_json(
        "fundamental_action_reconstruction/generated/p62_legacy_physical_role_transfer_to_strict_gate_kernel_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "f4_packet_present",
            "actual": f4["status"],
            "expected": "F4_EXECUTED_LEGACY_PHYSICAL_ROLE_TRANSFER_CLASSIFICATION_PACKET_NO_FALSE_PASS",
            "meaning": "F4 classifies the three legacy physical roles on the current repo state",
        },
        {
            "id": "p62_probe_negative",
            "actual": p62["status"],
            "expected": "CURRENT_REPO_EXPORTS_LEGACY_PHYSICAL_ROLE_CLAIMS_AND_STRICT_KERNEL_PIPELINE_BUT_NO_RIGOROUS_LEGACY_ROLE_TRANSFER_TO_STRICT_GATE_KERNEL_AFTER_P62",
            "meaning": "P62 confirms that no rigorous transfer of the three legacy physical roles is currently exported",
        },
        {
            "id": "weinberg_role_transfer_absent",
            "actual": p62["role_transfer_state"]["weinberg_angle_role_transfer_present"],
            "expected": False,
            "meaning": "the legacy Weinberg-angle role is not transferred onto K_strict_gate",
        },
        {
            "id": "fine_structure_role_transfer_absent",
            "actual": p62["role_transfer_state"]["fine_structure_role_transfer_present"],
            "expected": False,
            "meaning": "the legacy fine-structure role is not transferred onto K_strict_gate",
        },
        {
            "id": "gravity_hierarchy_role_transfer_absent",
            "actual": p62["role_transfer_state"]["gravity_hierarchy_role_transfer_present"],
            "expected": False,
            "meaning": "the legacy gravity-hierarchy role is not transferred onto K_strict_gate",
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
            "step": "N65",
            "status": "N65_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_PHYSICAL_ROLE_TRANSFER_STATE",
            "scope": "current_legacy_physical_role_transfer_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N65",
            "status": "N65_DISCHARGED_CURRENT_LEGACY_PHYSICAL_ROLE_TRANSFER_TO_STRICT_GATE_KERNEL_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_physical_role_transfer_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "legacy_weinberg_role_transfer_present": False,
                "legacy_fine_structure_role_transfer_present": False,
                "legacy_gravity_hierarchy_role_transfer_present": False,
                "rigorous_transfer_of_legacy_physical_roles_to_strict_gate_kernel_is_unsupported": True,
                "strict_operational_outputs_are_not_a_bridge_by_themselves": True,
                "full_closure_pass": False,
            },
            "missing_structure_classes": [
                "explicit_legacy_to_strict_kernel_bridge_operator_or_theorem",
                "explicit_claim_specific_role_transfer_theorem_for_legacy_weinberg_angle_semantics_onto_K_strict_gate",
                "explicit_claim_specific_role_transfer_theorem_for_legacy_fine_structure_semantics_onto_K_strict_gate",
                "explicit_claim_specific_role_transfer_theorem_for_legacy_gravity_hierarchy_semantics_onto_K_strict_gate",
                "explicit_retained_vs_replaced_partition_theorem_for_legacy_physical_roles_on_the_strict_side",
            ],
            "hard_limits": [
                "no_proof_that_the_strict_gate_kernel_is_false",
                "no_proof_that_the_legacy_formulas_are_false",
                "no_proof_that_no_future_role_transfer_can_ever_exist",
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
