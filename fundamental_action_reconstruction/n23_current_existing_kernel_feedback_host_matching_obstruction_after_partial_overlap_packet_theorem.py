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
        "R13": load_json("fundamental_action_reconstruction/generated/r13_partial_host_to_canonical_block_overlap_packet_summary.json"),
        "P19": load_json(
            "fundamental_action_reconstruction/generated/p19_existing_kernel_feedback_host_to_concrete_psi_block_identification_rerun_after_explicit_canonical_psi_block_export_summary.json"
        ),
        "P20": load_json(
            "fundamental_action_reconstruction/generated/p20_existing_kernel_feedback_host_matching_witness_rerun_after_partial_overlap_packet_summary.json"
        ),
        "QW2191": load_json(
            "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
        ),
        "C10": load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json"),
    }

    checks_spec = [
        {
            "id": "r13_partial_overlap_packet_present",
            "actual": sources["R13"]["status"],
            "expected": "PASS_PARTIAL_HOST_TO_CANONICAL_BLOCK_OVERLAP_PACKET_READY",
            "meaning": "R13 added the partial host/block overlap packet",
        },
        {
            "id": "p19_route_negative_before_r13",
            "actual": sources["P19"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_TO_CONCRETE_PSI_BLOCK_IDENTIFICATION_ROUTE_AFTER_R11_AND_R12_EXPLICIT_CANONICAL_PSI_BLOCK_EXPORT",
            "meaning": "before R13, the route still lacked a host-matching witness",
        },
        {
            "id": "p20_route_negative_after_r13",
            "actual": sources["P20"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_WITNESS_ROUTE_AFTER_R13_PARTIAL_OVERLAP_PACKET",
            "meaning": "after R13, the host-matching route still remains noncomputable",
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
        / "n23_current_existing_kernel_feedback_host_matching_obstruction_after_partial_overlap_packet_theorem_summary.json"
    )

    if mismatches:
        summary = {
            "step": "N23",
            "status": "N23_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_HOST_MATCHING_FRONTIER",
            "goal": "Check whether the current repo identifies the QW-2186 existing-feedback host with the exported canonical Psi block after adding the partial overlap packet.",
            "scope": "current_existing_kernel_feedback_host_matching_route_after_R13_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected host-matching frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_HOST_MATCHING_FRONTIER_BEFORE_CLAIMING_N23",
        }
    else:
        summary = {
            "step": "N23",
            "status": "N23_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_OBSTRUCTION_AFTER_R13_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after adding the partial host/canonical-block overlap packet, the current repo still does not identify the QW-2186 existing-feedback host with the exported canonical Psi block.",
            "scope": "current_existing_kernel_feedback_host_matching_route_after_R13_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "partial_host_to_canonical_block_overlap_present": True,
                "shared_kernel_light_facing_channel_present": True,
                "host_diagonal_floor_provenance_present": True,
                "full_physical_uniqueness_present": False,
                "host_matching_witness_present": False,
            },
            "missing_structure_classes": [
                "explicit_coefficient_specialization_witness_from_the_symbolic_canonical_kernel_channel_(K_i_j_plus_K_j_i)_over_2_to_the_frozen_numeric_K_total_matrix_on_the_same_12_slot_carrier",
                "explicit_diagonal_sector_equality_or_matching_witness_linking_the_host_floor_m0_squared_I_to_the_canonical_local_diagonal_sector_or_to_a_declared_control_pullback_of_it",
                "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future factorization routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "EITHER_PROVE_SELECTOR_RELEVANT_CANONICALIZATION_AND_ADD_KERNEL_SPECIALIZATION_PLUS_DIAGONAL_MATCHING_WITNESSES_OR_KEEP_THE_HOST_ROUTE_NEGATIVE",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
