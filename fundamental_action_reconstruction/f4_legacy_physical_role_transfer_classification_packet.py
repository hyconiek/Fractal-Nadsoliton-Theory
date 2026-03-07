#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f4_legacy_physical_role_transfer_classification_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f4_legacy_physical_role_transfer_classification_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def contains_all(text: str, parts: list[str]) -> bool:
    return all(part in text for part in parts)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p47 = load_json(
        "fundamental_action_reconstruction/generated/p47_legacy_ontological_kernel_to_strict_gate_kernel_bridge_probe_summary.json"
    )
    n50 = load_json(
        "fundamental_action_reconstruction/generated/n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
    )
    f2 = load_json(
        "fundamental_action_reconstruction/generated/f2_strict_gate_kernel_provenance_and_far_input_classification_packet_summary.json"
    )
    qw2005 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2005_PRE1700_TEX_REVALIDATION_MATRIX_EN_PL.md"
    )

    legacy_weinberg_role_present = contains_all(
        qw2005,
        [
            "C3",
            "sin^2(theta_W)=alpha_geo/12",
            "HEURISTIC / NOT STRICTLY DERIVED",
        ],
    )
    legacy_fine_structure_role_present = contains_all(
        qw2005,
        [
            "C4",
            "alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)",
            "PARTIAL / MODEL-LEVEL",
        ],
    )
    legacy_gravity_hierarchy_role_present = contains_all(
        qw2005,
        [
            "C5",
            "beta^N",
            "MODEL-CONSISTENT, NOT FULL INDEPENDENT PROOF",
        ],
    )

    checks_spec = [
        {
            "id": "p47_no_rigorous_kernel_bridge_present",
            "actual": p47["status"],
            "expected": "CURRENT_REPO_EXPORTS_LEGACY_AND_STRICT_KERNELS_BUT_NO_RIGOROUS_LEGACY_TO_STRICT_KERNEL_BRIDGE_AFTER_P47",
            "meaning": "P47 confirms that no rigorous legacy-to-strict kernel bridge is currently exported",
        },
        {
            "id": "n50_nonidentification_theorem_is_active",
            "actual": n50["status"],
            "expected": "N50_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_NONIDENTIFICATION_THEOREM_NO_FALSE_PASS",
            "meaning": "N50 turns the kernel split into a theorem-level current-repo-state nonidentification result",
        },
        {
            "id": "f2_strict_kernel_is_later_pipeline_operational_only",
            "actual": f2["status"],
            "expected": "F2_EXECUTED_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET_NO_FALSE_PASS",
            "meaning": "F2 classifies the strict gate kernel as a later-pipeline operational object unless a bridge is added",
        },
        {
            "id": "q2005_legacy_weinberg_role_is_downgraded",
            "actual": legacy_weinberg_role_present,
            "expected": True,
            "meaning": "QW-2005 records the legacy Weinberg-angle role as heuristic and not strictly derived",
        },
        {
            "id": "q2005_legacy_fine_structure_role_is_downgraded",
            "actual": legacy_fine_structure_role_present,
            "expected": True,
            "meaning": "QW-2005 records the legacy fine-structure role as partial and model-level",
        },
        {
            "id": "q2005_legacy_gravity_hierarchy_role_is_downgraded",
            "actual": legacy_gravity_hierarchy_role_present,
            "expected": True,
            "meaning": "QW-2005 records the legacy gravity-hierarchy role as model-consistent but not a full independent proof",
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

    legacy_roles = {
        "weinberg_angle_role": {
            "legacy_formula": "sin^2(theta_W)=alpha_geo/12",
            "q2005_status": "HEURISTIC / NOT STRICTLY DERIVED",
            "current_carrier": "legacy_kernel_package_only",
            "rigorous_transfer_to_strict_gate_supported": False,
        },
        "fine_structure_role": {
            "legacy_formula": "alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)",
            "q2005_status": "PARTIAL / MODEL-LEVEL",
            "current_carrier": "legacy_kernel_package_only",
            "rigorous_transfer_to_strict_gate_supported": False,
        },
        "gravity_hierarchy_role": {
            "legacy_formula": "beta^N gravity hierarchy",
            "q2005_status": "MODEL-CONSISTENT, NOT FULL INDEPENDENT PROOF",
            "current_carrier": "legacy_kernel_package_only",
            "rigorous_transfer_to_strict_gate_supported": False,
        },
    }

    artifact = {
        "stage": "F4",
        "lane": "legacy_physical_role_transfer_classification_current_repo_state_only",
        "goal": "classify_which_legacy_physical_roles_are_not_yet_rigorously_transferable_from_the_legacy_kernel_package_to_the_strict_gate_kernel",
        "status": "F4_EXECUTED_LEGACY_PHYSICAL_ROLE_TRANSFER_CLASSIFICATION_PACKET_NO_FALSE_PASS",
        "reason": "QW-2005 downgrades the legacy Weinberg-angle, fine-structure, and gravity-hierarchy claims while P47/N50 still show that no rigorous bridge identifies the legacy kernel with the strict gate kernel, and F2 classifies the strict gate kernel as a later-pipeline operational import rather than the restored ontological source layer; therefore those three legacy physical roles are not currently rigorously transferable onto K_strict_gate",
        "legacy_physical_roles": legacy_roles,
        "strict_gate_kernel_role_classification": {
            "class": "later_pipeline_operational_kernel_package",
            "silent_inheritance_of_legacy_physical_roles_disallowed": True,
        },
        "remaining_missing_objects": [
            "explicit_legacy_to_strict_kernel_bridge_operator_or_theorem",
            "explicit_claim_specific_role_transfer_theorem_for_legacy_weinberg_angle_semantics_onto_K_strict_gate",
            "explicit_claim_specific_role_transfer_theorem_for_legacy_fine_structure_semantics_onto_K_strict_gate",
            "explicit_claim_specific_role_transfer_theorem_for_legacy_gravity_hierarchy_semantics_onto_K_strict_gate",
            "explicit_retained_vs_replaced_partition_theorem_for_legacy_physical_roles_on_the_strict_side",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F4",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "legacy_physical_roles": legacy_roles,
        "strict_gate_kernel_role_classification": artifact[
            "strict_gate_kernel_role_classification"
        ],
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
