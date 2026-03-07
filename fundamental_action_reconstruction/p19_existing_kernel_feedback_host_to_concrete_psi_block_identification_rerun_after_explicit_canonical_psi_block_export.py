#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p19_existing_kernel_feedback_host_to_concrete_psi_block_identification_rerun_after_explicit_canonical_psi_block_export.json"
OUT_SUMMARY = GENERATED / "p19_existing_kernel_feedback_host_to_concrete_psi_block_identification_rerun_after_explicit_canonical_psi_block_export_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p18 = load_json(
        "fundamental_action_reconstruction/generated/p18_existing_kernel_feedback_host_to_concrete_psi_block_identification_rerun_after_symmetry_certified_transport_packet.json"
    )
    r11 = load_json(
        "fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json"
    )
    r12 = load_json(
        "fundamental_action_reconstruction/generated/r12_explicit_canonical_psi_block_export_packet_with_kernel_mixing_for_host_route.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")

    remaining_missing = [
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
        "explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_the_exported_canonical_Psi_block_on_full_declared_transport_support_or_its_declared_control_pullback",
    ]

    route_checks = [
        {
            "id": "p18_second_missing_object_was_concrete_submatrix_export",
            "pass": "explicit_assembled_and_coefficient_filled_concrete_Psi_sector_quadratic_submatrix_on_a_chosen_transported_index_set"
            in p18["remaining_missing_upstream_objects"],
            "expected": True,
            "actual": "explicit_assembled_and_coefficient_filled_concrete_Psi_sector_quadratic_submatrix_on_a_chosen_transported_index_set"
            in p18["remaining_missing_upstream_objects"],
            "meaning": "P18 still had the concrete Psi-sector submatrix export as a missing object",
        },
        {
            "id": "r11_explicit_declared_transport_packet_present",
            "pass": r11["result"]["explicit_declared_control_transport_packet_present"] is True,
            "expected": True,
            "actual": r11["result"]["explicit_declared_control_transport_packet_present"],
            "meaning": "the explicit declared transport packet is present",
        },
        {
            "id": "r12_explicit_canonical_block_export_present",
            "pass": r12["canonical_psi_block_export"]["shape"] == [12, 12],
            "expected": [12, 12],
            "actual": r12["canonical_psi_block_export"]["shape"],
            "meaning": "R12 exports a concrete coefficient-filled canonical Psi block",
        },
        {
            "id": "r12_light_facing_kernel_mixing_present",
            "pass": r12["kernel_mixing_light_facing_summary"]["off_diagonal_entry_count"] == 132,
            "expected": 132,
            "actual": r12["kernel_mixing_light_facing_summary"]["off_diagonal_entry_count"],
            "meaning": "the exported block contains the full off-diagonal kernel-mixing carrier",
        },
        {
            "id": "q2191_full_physical_uniqueness_still_open",
            "pass": q2191["flags"]["full_physical_uniqueness_closed"] is False,
            "expected": False,
            "actual": q2191["flags"]["full_physical_uniqueness_closed"],
            "meaning": "QW-2191 still blocks full physical uniqueness",
        },
        {
            "id": "c10_host_to_concrete_psi_block_identification_not_shown",
            "pass": c10["result"]["host_to_concrete_psi_block_identification"] == "not_shown",
            "expected": "not_shown",
            "actual": c10["result"]["host_to_concrete_psi_block_identification"],
            "meaning": "no host-to-concrete Psi-block identification is yet exported at audit level",
        },
    ]

    route_state = {
        "explicit_declared_control_transport_packet_present": True,
        "explicit_coefficient_filled_canonical_Psi_block_present": True,
        "light_facing_kernel_mixing_channel_present": True,
        "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        "host_to_submatrix_matching_witness_present": False,
    }

    report = {
        "stage": "P19",
        "goal": "rerun_compute_or_fail_existing_kernel_feedback_host_to_concrete_Psi_block_identification_after_R12",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_ROUTE_AFTER_R11_AND_R12_EXPLICIT_CANONICAL_PSI_BLOCK_EXPORT",
        "reason": "the repo now contains an explicit declared transport packet and a concrete coefficient-filled canonical Psi block with the full strict-admissible kernel-mixing carrier, but QW-2191 still leaves physical canonicalization of that transport open and the repo still exports no host-to-submatrix matching witness linking QW-2186 to the exported block or its declared control pullback",
        "lane": "existing_kernel_feedback_host_to_concrete_Psi_block_identification_route_after_R11_and_R12",
        "route_under_test": [
            "explicit_declared_control_transport_packet",
            "explicit_coefficient_filled_canonical_Psi_block_on_full_declared_transport_support",
            "light_facing_kernel_mixing_channel",
            "full_physical_uniqueness_or_selector_relevant_canonicalization",
            "host_to_submatrix_matching_witness",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R11_explicit_declared_control_transport_packet",
            "R12_explicit_canonical_Psi_block_export",
            "QW2166_exhaustive_canonical_hessian_gate",
        ],
        "decomposition_of_P18_second_missing_object": {
            "from_P18": "explicit_assembled_and_coefficient_filled_concrete_Psi_sector_quadratic_submatrix_on_a_chosen_transported_index_set",
            "discharged_by_R12_as": "explicit_coefficient_filled_canonical_Psi_block_on_full_declared_transport_support_with_kernel_mixing_highlight",
            "remaining_host_identification_gap": "explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_the_exported_canonical_Psi_block_on_full_declared_transport_support_or_its_declared_control_pullback",
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "blocking_frontier": {
            "QW2191_required_next_step": q2191["required_next_step"],
            "R11_B1": r11["frontier_text"],
            "R12_B1": r12["frontier_text"],
            "C10_B1": c10["residual_blockers"]["C10_B1"],
        },
        "required_next_step": "EITHER_PROVE_SELECTOR_RELEVANT_CANONICALIZATION_WITHIN_THE_QW2191_FAMILY_AND_FURNISH_A_HOST_TO_BLOCK_MATCHING_WITNESS_OR_KEEP_THE_HOST_IDENTIFICATION_ROUTE_NEGATIVE",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P19",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P18_second_missing_object": report["decomposition_of_P18_second_missing_object"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
