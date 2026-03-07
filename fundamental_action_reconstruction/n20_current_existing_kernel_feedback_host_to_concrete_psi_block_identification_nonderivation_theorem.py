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
        "P16": load_json("fundamental_action_reconstruction/generated/p16_existing_kernel_feedback_legacy_chart_reduced_operator_export_probe_summary.json"),
        "P17": load_json("fundamental_action_reconstruction/generated/p17_existing_kernel_feedback_host_to_concrete_psi_block_identification_probe_summary.json"),
        "C10": load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json"),
        "C11": load_json("fundamental_action_reconstruction/generated/c11_psi_sector_block_extraction_audit_summary.json"),
        "C12": load_json("fundamental_action_reconstruction/generated/c12_minimal_psi_block_extraction_packet_summary.json"),
        "C13": load_json("fundamental_action_reconstruction/generated/c13_mode_basis_control_index_set_audit_summary.json"),
        "C14": load_json("fundamental_action_reconstruction/generated/c14_control_mode_to_psi_transport_schema_summary.json"),
        "C20": load_json("fundamental_action_reconstruction/generated/c20_finite_materialization_recipe_audit_summary.json"),
    }

    checks_spec = [
        {
            "id": "p16_host_identification_was_the_first_missing_blocker",
            "actual": sources["P16"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_LEGACY_CHART_REDUCED_OPERATOR_EXPORT_ROUTE",
            "meaning": "P16 localized the first upstream blocker at host-to-concrete Psi-block identification",
        },
        {
            "id": "p17_host_to_concrete_psi_block_route_negative",
            "actual": sources["P17"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_ROUTE",
            "meaning": "the dedicated host-to-concrete Psi-block identification route remains noncomputable",
        },
        {
            "id": "c13_mode_basis_control_index_set_present",
            "actual": sources["C13"]["result"]["mode_basis_control_index_set_present"],
            "expected": "yes",
            "meaning": "a deterministic mode-basis control index-set is present",
        },
        {
            "id": "c14_control_transport_schema_present",
            "actual": sources["C14"]["result"]["control_transport_schema_present"],
            "expected": "yes",
            "meaning": "the control transport schema is present",
        },
        {
            "id": "c14_physical_transport_canonicalization_not_shown",
            "actual": sources["C14"]["result"]["strict_physical_justification_present"],
            "expected": "not_shown",
            "meaning": "the transport is not physically canonicalized",
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
            "meaning": "no assembled Psi-sector submatrix export is present",
        },
        {
            "id": "c20_serialization_run_not_executed",
            "actual": sources["C20"]["persisted_outputs"]["persisted_12_row_serialization_run_present"],
            "expected": False,
            "meaning": "the finite materialization recipe has not been executed into a persisted serialization run",
        },
        {
            "id": "c10_host_to_concrete_psi_block_identification_not_shown",
            "actual": sources["C10"]["result"]["host_to_concrete_psi_block_identification"],
            "expected": "not_shown",
            "meaning": "the host-to-concrete Psi-block identification remains absent",
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
        / "n20_current_existing_kernel_feedback_host_to_concrete_psi_block_identification_nonderivation_theorem_summary.json"
    )

    if mismatches:
        summary = {
            "step": "N20",
            "status": "N20_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_HOST_TO_CONCRETE_PSI_BLOCK_FRONTIER",
            "goal": "Check whether the current repo identifies the QW-2186 existing-feedback host operator with a concrete Psi-sector quadratic block.",
            "scope": "current_existing_kernel_feedback_host_to_concrete_psi_block_identification_route_after_P16_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected host-to-concrete-Psi-block frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_HOST_TO_CONCRETE_PSI_BLOCK_FRONTIER_BEFORE_CLAIMING_N20",
        }
    else:
        summary = {
            "step": "N20",
            "status": "N20_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_NONDERIVATION_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after deterministic control index-set declaration and control transport schema packetization, the current repo still does not identify the QW-2186 existing-feedback host operator with a concrete Psi-sector quadratic block.",
            "scope": "current_existing_kernel_feedback_host_to_concrete_psi_block_identification_route_after_P16_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "deterministic_mode_basis_control_index_set_present": True,
                "control_transport_schema_present": True,
                "concrete_psi_sector_block_extracted": False,
                "host_to_concrete_psi_block_identification_present": False,
            },
            "missing_structure_classes": [
                "strict_physical_canonicalization_of_the_control_transport_from_mode_basis_to_canonical_Psi_basis_for_selector_relevant_block_extraction",
                "explicit_assembled_and_coefficient_filled_concrete_Psi_sector_quadratic_submatrix_on_a_chosen_transported_index_set",
                "explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_that_concrete_Psi_sector_block",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future factorization routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "EITHER_EXPORT_A_PHYSICALLY_CANONICALIZED_CONCRETE_PSI_SECTOR_BLOCK_AND_MATCH_IT_TO_THE_QW2186_HOST_OR_KEEP_THE_HOST_IDENTIFICATION_ROUTE_NEGATIVE",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
