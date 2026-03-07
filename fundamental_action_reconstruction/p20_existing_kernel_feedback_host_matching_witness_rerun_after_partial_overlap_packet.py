#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p20_existing_kernel_feedback_host_matching_witness_rerun_after_partial_overlap_packet.json"
OUT_SUMMARY = GENERATED / "p20_existing_kernel_feedback_host_matching_witness_rerun_after_partial_overlap_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p19 = load_json(
        "fundamental_action_reconstruction/generated/p19_existing_kernel_feedback_host_to_concrete_psi_block_identification_rerun_after_explicit_canonical_psi_block_export.json"
    )
    r13 = load_json("fundamental_action_reconstruction/generated/r13_partial_host_to_canonical_block_overlap_packet.json")
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    remaining_missing = [
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
        "explicit_coefficient_specialization_witness_from_the_symbolic_canonical_kernel_channel_(K_i_j_plus_K_j_i)_over_2_to_the_frozen_numeric_K_total_matrix_on_the_same_12_slot_carrier",
        "explicit_diagonal_sector_equality_or_matching_witness_linking_the_host_floor_m0_squared_I_to_the_canonical_local_diagonal_sector_or_to_a_declared_control_pullback_of_it",
    ]

    route_checks = [
        {
            "id": "p19_matching_witness_was_still_missing",
            "pass": "explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_the_exported_canonical_Psi_block_on_full_declared_transport_support_or_its_declared_control_pullback"
            in p19["remaining_missing_upstream_objects"],
            "expected": True,
            "actual": "explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_the_exported_canonical_Psi_block_on_full_declared_transport_support_or_its_declared_control_pullback"
            in p19["remaining_missing_upstream_objects"],
            "meaning": "P19 still had the full host-matching witness as a missing object",
        },
        {
            "id": "r13_partial_overlap_present",
            "pass": r13["partial_overlap_witness"]["shared_12_slot_carrier_present"] is True,
            "expected": True,
            "actual": r13["partial_overlap_witness"]["shared_12_slot_carrier_present"],
            "meaning": "R13 adds a real partial host/block overlap packet",
        },
        {
            "id": "r13_shared_kernel_channel_present",
            "pass": r13["partial_overlap_witness"]["shared_off_diagonal_kernel_channel_present"] is True,
            "expected": True,
            "actual": r13["partial_overlap_witness"]["shared_off_diagonal_kernel_channel_present"],
            "meaning": "R13 makes the shared kernel/light-facing channel explicit",
        },
        {
            "id": "r13_host_diagonal_floor_provenance_present",
            "pass": r13["partial_overlap_witness"]["host_diagonal_floor_has_scalar_vacuum_provenance"] is True,
            "expected": True,
            "actual": r13["partial_overlap_witness"]["host_diagonal_floor_has_scalar_vacuum_provenance"],
            "meaning": "R13 makes the host diagonal floor provenance explicit",
        },
        {
            "id": "q2191_full_physical_uniqueness_still_open",
            "pass": q2191["flags"]["full_physical_uniqueness_closed"] is False,
            "expected": False,
            "actual": q2191["flags"]["full_physical_uniqueness_closed"],
            "meaning": "QW-2191 still blocks full physical uniqueness",
        },
        {
            "id": "c10_host_identification_not_shown",
            "pass": c10["result"]["host_to_concrete_psi_block_identification"] == "not_shown",
            "expected": "not_shown",
            "actual": c10["result"]["host_to_concrete_psi_block_identification"],
            "meaning": "host-to-concrete Psi-block identification remains absent at audit level",
        },
    ]

    route_state = {
        "partial_host_to_canonical_block_overlap_present": True,
        "shared_kernel_light_facing_channel_present": True,
        "host_diagonal_floor_provenance_present": True,
        "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        "kernel_coefficient_specialization_witness_present": False,
        "diagonal_sector_matching_witness_present": False,
    }

    report = {
        "stage": "P20",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_matching_witness_after_R13",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_WITNESS_ROUTE_AFTER_R13_PARTIAL_OVERLAP_PACKET",
        "reason": "the repo now contains a real partial overlap packet between the QW-2186 host operator and the exported canonical Psi block, including the shared kernel/light-facing channel and host diagonal-floor provenance, but it still lacks the coefficient-specialization witness for the kernel part, the diagonal-sector matching witness, and the QW-2191 canonicalization boundary remains open",
        "lane": "existing_kernel_feedback_host_matching_witness_route_after_R13",
        "route_under_test": [
            "partial_host_to_canonical_block_overlap_packet",
            "kernel_coefficient_specialization_witness",
            "diagonal_sector_matching_witness",
            "full_physical_uniqueness_or_selector_relevant_canonicalization",
            "host_to_concrete_Psi_block_identification",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R13_partial_host_to_canonical_block_overlap_packet",
            "R12_explicit_canonical_Psi_block_export",
            "R8_existing_kernel_feedback_host_carrier_packet",
        ],
        "decomposition_of_P19_matching_gap": {
            "from_P19": "explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_the_exported_canonical_Psi_block_on_full_declared_transport_support_or_its_declared_control_pullback",
            "into_current_blockers": remaining_missing[1:],
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "QW2191_required_next_step": q2191["required_next_step"],
            "R13_B1": r13["frontier_text"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_SELECTOR_RELEVANT_CANONICALIZATION_AND_ADD_KERNEL_SPECIALIZATION_PLUS_DIAGONAL_MATCHING_WITNESSES_OR_KEEP_THE_HOST_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P20",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P19_matching_gap": report["decomposition_of_P19_matching_gap"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
