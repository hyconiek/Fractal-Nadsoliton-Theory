#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p17_existing_kernel_feedback_host_to_concrete_psi_block_identification_probe.json"
OUT_SUMMARY = GENERATED / "p17_existing_kernel_feedback_host_to_concrete_psi_block_identification_probe_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p16 = load_json(
        "fundamental_action_reconstruction/generated/p16_existing_kernel_feedback_legacy_chart_reduced_operator_export_probe.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")
    c11 = load_json("fundamental_action_reconstruction/generated/c11_psi_sector_block_extraction_audit_summary.json")
    c12 = load_json("fundamental_action_reconstruction/generated/c12_minimal_psi_block_extraction_packet_summary.json")
    c13 = load_json("fundamental_action_reconstruction/generated/c13_mode_basis_control_index_set_audit_summary.json")
    c14 = load_json("fundamental_action_reconstruction/generated/c14_control_mode_to_psi_transport_schema_summary.json")
    c20 = load_json("fundamental_action_reconstruction/generated/c20_finite_materialization_recipe_audit_summary.json")

    p16_missing = p16["remaining_missing_upstream_objects"]
    remaining_missing = [
        "strict_physical_canonicalization_of_the_control_transport_from_mode_basis_to_canonical_Psi_basis_for_selector_relevant_block_extraction",
        "explicit_assembled_and_coefficient_filled_concrete_Psi_sector_quadratic_submatrix_on_a_chosen_transported_index_set",
        "explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_that_concrete_Psi_sector_block",
    ]

    route_checks = [
        {
            "id": "p16_originally_missing_host_to_concrete_psi_block_identification",
            "pass": "host_to_concrete_Psi_sector_quadratic_block_identification_for_the_existing_kernel_feedback_host_operator"
            in p16_missing,
            "expected": True,
            "actual": "host_to_concrete_Psi_sector_quadratic_block_identification_for_the_existing_kernel_feedback_host_operator"
            in p16_missing,
            "meaning": "P16 indeed localized the first upstream blocker as host-to-concrete Psi-block identification",
        },
        {
            "id": "c13_mode_basis_control_index_set_present",
            "pass": c13["result"]["mode_basis_control_index_set_present"] == "yes",
            "expected": "yes",
            "actual": c13["result"]["mode_basis_control_index_set_present"],
            "meaning": "the deterministic control index-set exists in mode basis",
        },
        {
            "id": "c14_control_transport_schema_present",
            "pass": c14["result"]["control_transport_schema_present"] == "yes",
            "expected": "yes",
            "actual": c14["result"]["control_transport_schema_present"],
            "meaning": "a control transport schema from mode basis to canonical Psi carrier is present",
        },
        {
            "id": "c14_strict_physical_justification_not_shown",
            "pass": c14["result"]["strict_physical_justification_present"] == "not_shown",
            "expected": "not_shown",
            "actual": c14["result"]["strict_physical_justification_present"],
            "meaning": "the transport has not been physically canonicalized for selector-relevant extraction",
        },
        {
            "id": "c11_concrete_psi_block_not_extracted",
            "pass": c11["result"]["concrete_block_extracted"] == "not_shown",
            "expected": "not_shown",
            "actual": c11["result"]["concrete_block_extracted"],
            "meaning": "no concrete Psi-sector block has been explicitly extracted",
        },
        {
            "id": "c12_assembled_submatrix_not_exported",
            "pass": c12["result"]["assembled_submatrix_exported"] == "not_shown",
            "expected": "not_shown",
            "actual": c12["result"]["assembled_submatrix_exported"],
            "meaning": "no assembled Psi-sector submatrix export is present",
        },
        {
            "id": "c20_executed_serialization_run_not_shown",
            "pass": c20["persisted_outputs"]["persisted_12_row_serialization_run_present"] is False,
            "expected": False,
            "actual": c20["persisted_outputs"]["persisted_12_row_serialization_run_present"],
            "meaning": "the finite materialization recipe still has no executed persisted serialization run",
        },
        {
            "id": "c10_host_to_concrete_psi_block_identification_not_shown",
            "pass": c10["result"]["host_to_concrete_psi_block_identification"] == "not_shown",
            "expected": "not_shown",
            "actual": c10["result"]["host_to_concrete_psi_block_identification"],
            "meaning": "the host-to-concrete Psi-block identification is still absent at audit level",
        },
    ]

    route_state = {
        "deterministic_mode_basis_control_index_set_present": True,
        "control_transport_schema_present": True,
        "strict_physical_transport_canonicalization_present": False,
        "concrete_psi_sector_block_extracted": False,
        "coefficient_filled_concrete_psi_sector_submatrix_present": False,
        "host_to_submatrix_matching_witness_present": False,
    }

    report = {
        "stage": "P17",
        "goal": "compute_or_fail_existing_kernel_feedback_host_to_concrete_Psi_sector_block_identification",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_ROUTE",
        "reason": "the repo now contains a deterministic control index-set in mode basis and a control-transport schema into the canonical Psi carrier, but it still exports neither a physically canonicalized transport for selector-relevant block extraction, nor an assembled coefficient-filled concrete Psi-sector submatrix on a chosen transported index-set, nor a host-to-submatrix matching witness identifying the QW-2186 host operator with that concrete block",
        "lane": "existing_kernel_feedback_host_to_concrete_Psi_block_identification_route_after_P16",
        "route_under_test": [
            "deterministic_mode_basis_control_index_set",
            "control_transport_schema_to_canonical_Psi_basis",
            "strict_physical_transport_canonicalization",
            "concrete_Psi_sector_block_extraction",
            "coefficient_filled_concrete_Psi_sector_submatrix_export",
            "host_to_submatrix_matching_witness",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "C13_mode_basis_control_index_sets",
            "C14_control_transport_schema",
            "C12_minimal_Psi_block_extraction_packet",
            "C20_finite_materialization_recipe",
        ],
        "decomposition_of_P16_missing_object": {
            "from_P16": "host_to_concrete_Psi_sector_quadratic_block_identification_for_the_existing_kernel_feedback_host_operator",
            "into_current_blockers": remaining_missing,
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "C14_B1": c14["residual_blockers"]["C14_B1"],
            "C14_B2": c14["residual_blockers"]["C14_B2"],
            "C20_B1": c20["residual_blockers"]["C20_B1"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_EXPORT_A_PHYSICALLY_CANONICALIZED_CONCRETE_PSI_SECTOR_BLOCK_AND_MATCH_IT_TO_THE_QW2186_HOST_OR_KEEP_THE_HOST_IDENTIFICATION_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P17",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P16_missing_object": report["decomposition_of_P16_missing_object"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
