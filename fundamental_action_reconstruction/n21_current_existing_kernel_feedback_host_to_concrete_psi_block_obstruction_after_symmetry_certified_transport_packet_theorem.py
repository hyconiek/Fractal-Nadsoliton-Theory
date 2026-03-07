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
        "P17": load_json("fundamental_action_reconstruction/generated/p17_existing_kernel_feedback_host_to_concrete_psi_block_identification_probe_summary.json"),
        "R11": load_json("fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route_summary.json"),
        "P18": load_json(
            "fundamental_action_reconstruction/generated/p18_existing_kernel_feedback_host_to_concrete_psi_block_identification_rerun_after_symmetry_certified_transport_packet_summary.json"
        ),
        "QW2191": load_json(
            "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
        ),
        "C10": load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json"),
        "C11": load_json("fundamental_action_reconstruction/generated/c11_psi_sector_block_extraction_audit_summary.json"),
        "C12": load_json("fundamental_action_reconstruction/generated/c12_minimal_psi_block_extraction_packet_summary.json"),
        "C20": load_json("fundamental_action_reconstruction/generated/c20_finite_materialization_recipe_audit_summary.json"),
    }

    checks_spec = [
        {
            "id": "p17_route_negative_before_r11",
            "actual": sources["P17"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_ROUTE",
            "meaning": "P17 localized the route before the symmetry-certified transport packet was added",
        },
        {
            "id": "r11_partial_transport_packet_present",
            "actual": sources["R11"]["status"],
            "expected": "PASS_PARTIAL_EXPLICIT_DECLARED_CONTROL_TRANSPORT_PACKET_READY_BUT_PHYSICAL_UNIQUENESS_OPEN",
            "meaning": "R11 added the explicit declared transport packet with symmetry certification but without full physical uniqueness",
        },
        {
            "id": "p18_route_negative_after_r11",
            "actual": sources["P18"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_ROUTE_AFTER_R11_SYMMETRY_CERTIFIED_TRANSPORT_PACKET",
            "meaning": "the host-identification route still remains noncomputable after R11",
        },
        {
            "id": "qw2191_full_physical_uniqueness_closed_false",
            "actual": sources["QW2191"]["flags"]["full_physical_uniqueness_closed"],
            "expected": False,
            "meaning": "QW-2191 still blocks full physical uniqueness",
        },
        {
            "id": "c11_concrete_block_not_extracted",
            "actual": sources["C11"]["result"]["concrete_block_extracted"],
            "expected": "not_shown",
            "meaning": "no concrete Psi-sector block has been extracted",
        },
        {
            "id": "c12_assembled_submatrix_not_exported",
            "actual": sources["C12"]["result"]["assembled_submatrix_exported"],
            "expected": "not_shown",
            "meaning": "no coefficient-filled concrete Psi-sector submatrix export is present",
        },
        {
            "id": "c20_serialization_run_not_executed",
            "actual": sources["C20"]["persisted_outputs"]["persisted_12_row_serialization_run_present"],
            "expected": False,
            "meaning": "the finite materialization recipe is still not executed",
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
        / "n21_current_existing_kernel_feedback_host_to_concrete_psi_block_obstruction_after_symmetry_certified_transport_packet_theorem_summary.json"
    )

    if mismatches:
        summary = {
            "step": "N21",
            "status": "N21_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_TRANSPORT_CANONICALIZATION_FRONTIER",
            "goal": "Check whether the current repo identifies the QW-2186 existing-feedback host with a concrete Psi-sector quadratic block after adding the symmetry-certified declared transport packet.",
            "scope": "current_existing_kernel_feedback_host_to_concrete_psi_block_route_after_R11_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected transport-canonicalization frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_TRANSPORT_CANONICALIZATION_FRONTIER_BEFORE_CLAIMING_N21",
        }
    else:
        summary = {
            "step": "N21",
            "status": "N21_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_OBSTRUCTION_AFTER_R11_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after adding the explicit declared control transport packet and symmetry certification, the current repo still does not identify the QW-2186 existing-feedback host operator with a concrete Psi-sector quadratic block.",
            "scope": "current_existing_kernel_feedback_host_to_concrete_psi_block_route_after_R11_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "explicit_declared_transport_packet_present": True,
                "symmetry_certificate_present": True,
                "full_physical_uniqueness_present": False,
                "concrete_psi_sector_block_extracted": False,
                "host_to_concrete_psi_block_identification_present": False,
            },
            "missing_structure_classes": [
                "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
                "explicit_assembled_and_coefficient_filled_concrete_Psi_sector_quadratic_submatrix_on_a_chosen_transported_index_set",
                "explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_that_concrete_block",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future factorization routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "EITHER_PROVE_SELECTOR_RELEVANT_CANONICALIZATION_WITHIN_THE_QW2191_FAMILY_AND_EXPORT_A_CONCRETE_PSI_BLOCK_WITH_A_HOST_MATCH_OR_KEEP_THE_HOST_IDENTIFICATION_ROUTE_NEGATIVE",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
