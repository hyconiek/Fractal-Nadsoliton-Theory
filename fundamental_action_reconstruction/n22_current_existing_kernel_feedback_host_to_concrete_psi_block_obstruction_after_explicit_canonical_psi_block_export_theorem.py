#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "R11": load_json("fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route_summary.json"),
        "R12": load_json(
            "fundamental_action_reconstruction/generated/r12_explicit_canonical_psi_block_export_packet_with_kernel_mixing_for_host_route_summary.json"
        ),
        "P18": load_json(
            "fundamental_action_reconstruction/generated/p18_existing_kernel_feedback_host_to_concrete_psi_block_identification_rerun_after_symmetry_certified_transport_packet_summary.json"
        ),
        "P19": load_json(
            "fundamental_action_reconstruction/generated/p19_existing_kernel_feedback_host_to_concrete_psi_block_identification_rerun_after_explicit_canonical_psi_block_export_summary.json"
        ),
        "QW2191": load_json(
            "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
        ),
        "C10": load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json"),
    }

    checks_spec = [
        {
            "id": "r11_transport_packet_present",
            "actual": sources["R11"]["status"],
            "expected": "PASS_PARTIAL_EXPLICIT_DECLARED_CONTROL_TRANSPORT_PACKET_READY_BUT_PHYSICAL_UNIQUENESS_OPEN",
            "meaning": "R11 added the explicit declared control transport packet",
        },
        {
            "id": "p18_route_negative_after_r11",
            "actual": sources["P18"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_ROUTE_AFTER_R11_SYMMETRY_CERTIFIED_TRANSPORT_PACKET",
            "meaning": "the host-identification route remained negative after R11",
        },
        {
            "id": "r12_canonical_block_export_present",
            "actual": sources["R12"]["status"],
            "expected": "PASS_PARTIAL_EXPLICIT_CANONICAL_PSI_BLOCK_EXPORT_READY_ON_FULL_DECLARED_TRANSPORT_SUPPORT",
            "meaning": "R12 added the explicit canonical Psi block export",
        },
        {
            "id": "p19_route_negative_after_r12",
            "actual": sources["P19"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_ROUTE_AFTER_R11_AND_R12_EXPLICIT_CANONICAL_PSI_BLOCK_EXPORT",
            "meaning": "the route still remains negative after R12",
        },
        {
            "id": "qw2191_full_physical_uniqueness_closed_false",
            "actual": sources["QW2191"]["flags"]["full_physical_uniqueness_closed"],
            "expected": False,
            "meaning": "QW-2191 still blocks full physical uniqueness",
        },
        {
            "id": "c10_host_identification_not_shown",
            "actual": sources["C10"]["result"]["host_to_concrete_psi_block_identification"],
            "expected": "not_shown",
            "meaning": "host-to-concrete Psi-block identification remains absent",
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

    out = (
        ROOT
        / "generated"
        / "n22_current_existing_kernel_feedback_host_to_concrete_psi_block_obstruction_after_explicit_canonical_psi_block_export_theorem_summary.json"
    )

    if mismatches:
        summary = {
            "step": "N22",
            "status": "N22_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CANONICAL_BLOCK_EXPORT_FRONTIER",
            "goal": "Check whether the current repo identifies the QW-2186 existing-feedback host with a concrete Psi-sector quadratic block after adding the explicit canonical Psi block export.",
            "scope": "current_existing_kernel_feedback_host_to_concrete_psi_block_route_after_R11_and_R12_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected canonical-block-export frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_CANONICAL_BLOCK_EXPORT_FRONTIER_BEFORE_CLAIMING_N22",
        }
    else:
        summary = {
            "step": "N22",
            "status": "N22_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_OBSTRUCTION_AFTER_R12_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after exporting an explicit canonical Psi block with the strict-admissible kernel-mixing carrier, the current repo still does not identify the QW-2186 existing-feedback host operator with a concrete Psi-sector quadratic block.",
            "scope": "current_existing_kernel_feedback_host_to_concrete_psi_block_route_after_R11_and_R12_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "explicit_declared_transport_packet_present": True,
                "explicit_canonical_psi_block_export_present": True,
                "strict_admissible_light_facing_kernel_mixing_present": True,
                "full_physical_uniqueness_present": False,
                "host_to_concrete_psi_block_identification_present": False,
            },
            "missing_structure_classes": [
                "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
                "explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_the_exported_canonical_Psi_block_on_full_declared_transport_support_or_its_declared_control_pullback",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future factorization routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "EITHER_PROVE_SELECTOR_RELEVANT_CANONICALIZATION_WITHIN_THE_QW2191_FAMILY_AND_FURNISH_A_HOST_TO_BLOCK_MATCHING_WITNESS_OR_KEEP_THE_HOST_IDENTIFICATION_ROUTE_NEGATIVE",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
