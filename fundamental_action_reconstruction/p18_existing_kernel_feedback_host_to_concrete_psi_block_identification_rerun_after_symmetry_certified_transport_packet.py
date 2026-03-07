#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p18_existing_kernel_feedback_host_to_concrete_psi_block_identification_rerun_after_symmetry_certified_transport_packet.json"
OUT_SUMMARY = GENERATED / "p18_existing_kernel_feedback_host_to_concrete_psi_block_identification_rerun_after_symmetry_certified_transport_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p17 = load_json(
        "fundamental_action_reconstruction/generated/p17_existing_kernel_feedback_host_to_concrete_psi_block_identification_probe.json"
    )
    r11 = load_json(
        "fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")
    c11 = load_json("fundamental_action_reconstruction/generated/c11_psi_sector_block_extraction_audit_summary.json")
    c12 = load_json("fundamental_action_reconstruction/generated/c12_minimal_psi_block_extraction_packet_summary.json")
    c20 = load_json("fundamental_action_reconstruction/generated/c20_finite_materialization_recipe_audit_summary.json")

    remaining_missing = [
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
        "explicit_assembled_and_coefficient_filled_concrete_Psi_sector_quadratic_submatrix_on_a_chosen_transported_index_set",
        "explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_that_concrete_block",
    ]

    route_checks = [
        {
            "id": "p17_originally_missing_transport_canonicalization",
            "pass": "strict_physical_canonicalization_of_the_control_transport_from_mode_basis_to_canonical_Psi_basis_for_selector_relevant_block_extraction"
            in p17["remaining_missing_upstream_objects"],
            "expected": True,
            "actual": "strict_physical_canonicalization_of_the_control_transport_from_mode_basis_to_canonical_Psi_basis_for_selector_relevant_block_extraction"
            in p17["remaining_missing_upstream_objects"],
            "meaning": "P17 indeed localized physical transport canonicalization as the first missing structure class",
        },
        {
            "id": "r11_explicit_declared_transport_packet_present",
            "pass": r11["result"]["explicit_declared_control_transport_packet_present"] is True,
            "expected": True,
            "actual": r11["result"]["explicit_declared_control_transport_packet_present"],
            "meaning": "R11 adds the explicit declared control transport packet",
        },
        {
            "id": "r11_symmetry_certificate_present",
            "pass": r11["result"]["symmetry_certificate_present"] is True,
            "expected": True,
            "actual": r11["result"]["symmetry_certificate_present"],
            "meaning": "R11 adds the symmetry certificate for the declared transport",
        },
        {
            "id": "q2191_full_physical_uniqueness_still_open",
            "pass": q2191["flags"]["full_physical_uniqueness_closed"] is False,
            "expected": False,
            "actual": q2191["flags"]["full_physical_uniqueness_closed"],
            "meaning": "QW-2191 still blocks full physical uniqueness of mode-index assignment",
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
            "meaning": "no assembled coefficient-filled Psi-sector submatrix export is present",
        },
        {
            "id": "c20_serialization_run_not_executed",
            "pass": c20["persisted_outputs"]["persisted_12_row_serialization_run_present"] is False,
            "expected": False,
            "actual": c20["persisted_outputs"]["persisted_12_row_serialization_run_present"],
            "meaning": "the finite materialization recipe is still not executed into a persisted serialization run",
        },
        {
            "id": "c10_host_to_concrete_psi_block_identification_not_shown",
            "pass": c10["result"]["host_to_concrete_psi_block_identification"] == "not_shown",
            "expected": "not_shown",
            "actual": c10["result"]["host_to_concrete_psi_block_identification"],
            "meaning": "the host-to-concrete Psi-block identification remains absent",
        },
    ]

    route_state = {
        "explicit_declared_control_transport_packet_present": True,
        "symmetry_certificate_present": True,
        "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        "concrete_psi_sector_block_extracted": False,
        "coefficient_filled_concrete_psi_sector_submatrix_present": False,
        "host_to_submatrix_matching_witness_present": False,
    }

    report = {
        "stage": "P18",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_to_concrete_Psi_sector_block_identification_after_R11",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_ROUTE_AFTER_R11_SYMMETRY_CERTIFIED_TRANSPORT_PACKET",
        "reason": "the repo now contains a deterministic control index-set, a schema-level transport, an explicit declared transport matrix, and a symmetry certificate for that transport, but QW-2191 still leaves full physical uniqueness and selector-relevant canonicalization open, and the repo still exports neither a concrete coefficient-filled Psi-sector submatrix on a chosen transported index-set nor a host-to-submatrix matching witness",
        "lane": "existing_kernel_feedback_host_to_concrete_Psi_block_identification_route_after_R11",
        "route_under_test": [
            "deterministic_mode_basis_control_index_set",
            "control_transport_schema_to_canonical_Psi_basis",
            "explicit_declared_control_transport_packet",
            "symmetry_certificate_for_the_declared_transport",
            "full_physical_uniqueness_or_selector_relevant_canonicalization",
            "concrete_Psi_sector_block_extraction",
            "coefficient_filled_concrete_Psi_sector_submatrix_export",
            "host_to_submatrix_matching_witness",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "C13_mode_basis_control_index_sets",
            "C14_control_transport_schema",
            "R11_explicit_declared_control_transport_packet",
            "QW2190_symmetry_certificate",
        ],
        "decomposition_of_P17_first_missing_object": {
            "from_P17": "strict_physical_canonicalization_of_the_control_transport_from_mode_basis_to_canonical_Psi_basis_for_selector_relevant_block_extraction",
            "present_partial_component": "explicit_declared_control_transport_packet_and_symmetry_certificate_on_the_canonical_Psi_carrier",
            "remaining_missing_component": "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "R11_B1": r11["frontier_text"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
            "C14_B2": "no_assembled_Psi_x_Psi_submatrix_after_adopting_the_control_transport_schema",
            "QW2191_required_next_step": q2191["required_next_step"],
        },
        "required_next_step": "EITHER_PROVE_SELECTOR_RELEVANT_CANONICALIZATION_WITHIN_THE_QW2191_FAMILY_AND_EXPORT_A_CONCRETE_PSI_BLOCK_WITH_A_HOST_MATCH_OR_KEEP_THE_HOST_IDENTIFICATION_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P18",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P17_first_missing_object": report["decomposition_of_P17_first_missing_object"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
