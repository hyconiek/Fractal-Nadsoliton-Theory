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
    / "p62_legacy_physical_role_transfer_to_strict_gate_kernel_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p62_legacy_physical_role_transfer_to_strict_gate_kernel_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f4 = load_json(
        "fundamental_action_reconstruction/generated/f4_legacy_physical_role_transfer_classification_packet_summary.json"
    )
    p47 = load_json(
        "fundamental_action_reconstruction/generated/p47_legacy_ontological_kernel_to_strict_gate_kernel_bridge_probe_summary.json"
    )
    n50 = load_json(
        "fundamental_action_reconstruction/generated/n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
    )
    f2 = load_json(
        "fundamental_action_reconstruction/generated/f2_strict_gate_kernel_provenance_and_far_input_classification_packet_summary.json"
    )

    legacy_roles = f4["legacy_physical_roles"]

    checks_spec = [
        {
            "id": "f4_role_transfer_classification_packet_present",
            "actual": f4["status"],
            "expected": "F4_EXECUTED_LEGACY_PHYSICAL_ROLE_TRANSFER_CLASSIFICATION_PACKET_NO_FALSE_PASS",
            "meaning": "F4 classifies the legacy physical roles and their non-transfer status on the current repo state",
        },
        {
            "id": "p47_kernel_bridge_absent",
            "actual": p47["status"],
            "expected": "CURRENT_REPO_EXPORTS_LEGACY_AND_STRICT_KERNELS_BUT_NO_RIGOROUS_LEGACY_TO_STRICT_KERNEL_BRIDGE_AFTER_P47",
            "meaning": "P47 confirms that the kernel-level bridge is still absent",
        },
        {
            "id": "n50_nonidentification_theorem_active",
            "actual": n50["status"],
            "expected": "N50_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_NONIDENTIFICATION_THEOREM_NO_FALSE_PASS",
            "meaning": "N50 confirms that no rigorous legacy-to-strict kernel identification exists on the current repo state",
        },
        {
            "id": "f2_strict_kernel_still_operational_only",
            "actual": f2["status"],
            "expected": "F2_EXECUTED_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET_NO_FALSE_PASS",
            "meaning": "F2 still classifies the strict gate kernel as a later-pipeline operational object",
        },
        {
            "id": "weinberg_role_transfer_supported",
            "actual": legacy_roles["weinberg_angle_role"][
                "rigorous_transfer_to_strict_gate_supported"
            ],
            "expected": False,
            "meaning": "the legacy Weinberg-angle role is not rigorously transferred onto K_strict_gate",
        },
        {
            "id": "fine_structure_role_transfer_supported",
            "actual": legacy_roles["fine_structure_role"][
                "rigorous_transfer_to_strict_gate_supported"
            ],
            "expected": False,
            "meaning": "the legacy fine-structure role is not rigorously transferred onto K_strict_gate",
        },
        {
            "id": "gravity_hierarchy_role_transfer_supported",
            "actual": legacy_roles["gravity_hierarchy_role"][
                "rigorous_transfer_to_strict_gate_supported"
            ],
            "expected": False,
            "meaning": "the legacy gravity-hierarchy role is not rigorously transferred onto K_strict_gate",
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
        "explicit_legacy_to_strict_kernel_bridge_operator_or_theorem",
        "explicit_claim_specific_role_transfer_theorem_for_legacy_weinberg_angle_semantics_onto_K_strict_gate",
        "explicit_claim_specific_role_transfer_theorem_for_legacy_fine_structure_semantics_onto_K_strict_gate",
        "explicit_claim_specific_role_transfer_theorem_for_legacy_gravity_hierarchy_semantics_onto_K_strict_gate",
        "explicit_retained_vs_replaced_partition_theorem_for_legacy_physical_roles_on_the_strict_side",
    ]

    artifact = {
        "stage": "P62",
        "lane": "legacy_physical_role_transfer_to_strict_gate_kernel_probe_current_repo_state_only",
        "goal": "test_whether_the_current_repo_already_exports_a_rigorous_transfer_of_legacy_weinberg_fine_structure_and_gravity_hierarchy_roles_onto_the_strict_gate_kernel",
        "status": "CURRENT_REPO_EXPORTS_LEGACY_PHYSICAL_ROLE_CLAIMS_AND_STRICT_KERNEL_PIPELINE_BUT_NO_RIGOROUS_LEGACY_ROLE_TRANSFER_TO_STRICT_GATE_KERNEL_AFTER_P62",
        "reason": "the current repo exports the legacy physical-role package and the later strict kernel pipeline, but F4 keeps the three legacy physical roles on the legacy side only, P47/N50 keep the kernel bridge absent, and F2 keeps the strict gate kernel operational-only; therefore no rigorous role transfer from the legacy package to K_strict_gate is currently exported",
        "legacy_roles_checked": legacy_roles,
        "role_transfer_state": {
            "weinberg_angle_role_transfer_present": False,
            "fine_structure_role_transfer_present": False,
            "gravity_hierarchy_role_transfer_present": False,
        },
        "remaining_missing_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P62",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "role_transfer_state": artifact["role_transfer_state"],
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
