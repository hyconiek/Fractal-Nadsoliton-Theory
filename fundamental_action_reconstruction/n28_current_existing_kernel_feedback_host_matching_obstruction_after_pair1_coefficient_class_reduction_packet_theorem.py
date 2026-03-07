#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n28_current_existing_kernel_feedback_host_matching_obstruction_after_pair1_coefficient_class_reduction_packet_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "R18": load_json(
            "fundamental_action_reconstruction/generated/r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route_summary.json"
        ),
        "P24": load_json(
            "fundamental_action_reconstruction/generated/p24_existing_kernel_feedback_host_matching_witness_rerun_after_host_side_residual_absence_packet_summary.json"
        ),
        "P25": load_json(
            "fundamental_action_reconstruction/generated/p25_existing_kernel_feedback_host_matching_witness_rerun_after_pair1_coefficient_class_reduction_packet_summary.json"
        ),
        "QW2191": load_json(
            "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
        ),
        "C10": load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json"),
    }

    checks_spec = [
        {
            "id": "r18_pair1_reduction_packet_present",
            "actual": sources["R18"]["status"],
            "expected": "PASS_PARTIAL_PAIR1_RESIDUAL_DECLARED_PULLBACK_COEFFICIENT_CLASS_REDUCTION_PACKET_READY",
            "meaning": "R18 added the exact pair1 coefficient-class reduction packet",
        },
        {
            "id": "p24_route_negative_before_r18",
            "actual": sources["P24"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_WITNESS_ROUTE_AFTER_R17_HOST_SIDE_RESIDUAL_ABSENCE_PACKET",
            "meaning": "before R18, the route still lacked a generic zero witness for the declared residual pullback",
        },
        {
            "id": "p25_route_negative_after_r18",
            "actual": sources["P25"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_WITNESS_ROUTE_AFTER_R18_PAIR1_COEFFICIENT_CLASS_REDUCTION_PACKET",
            "meaning": "after R18, the host-matching route still remains noncomputable",
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

    if mismatches:
        summary = {
            "step": "N28",
            "status": "N28_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_HOST_MATCHING_FRONTIER",
            "goal": "Check whether the current repo identifies the QW-2186 existing-feedback host with the exported canonical Psi block after adding the pair1 coefficient-class reduction packet.",
            "scope": "current_existing_kernel_feedback_host_matching_route_after_R18_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected host-matching frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_HOST_MATCHING_FRONTIER_BEFORE_CLAIMING_N28",
        }
    else:
        summary = {
            "step": "N28",
            "status": "N28_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_OBSTRUCTION_AFTER_R18_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after adding the exact pair1 coefficient-class reduction packet, the current repo still does not identify the QW-2186 existing-feedback host with the exported canonical Psi block.",
            "scope": "current_existing_kernel_feedback_host_matching_route_after_R18_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "shared_kernel_light_channel_specialized": True,
                "host_scalar_floor_embedding_present": True,
                "declared_control_pullback_of_residual_local_diagonal_sector_present": True,
                "host_side_residual_diagonal_correction_present": False,
                "pair1_coefficient_class_reduction_present": True,
                "pair1_c1c1_zero_witness_present": False,
                "pair1_c1s1_zero_witness_present": False,
                "pair1_s1s1_zero_witness_present": False,
                "full_physical_uniqueness_present": False,
                "host_matching_witness_present": False,
            },
            "missing_structure_classes": [
                "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
                "explicit_zero_witness_for_the_declared_pair1_residual_c1s1_equation",
                "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
                "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future factorization routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "EITHER_PROVE_THE_THREE_INDEPENDENT_DECLARED_PAIR1_ZERO_EQUATIONS_AND_DISCHARGE_SELECTOR_RELEVANT_CANONICALIZATION_OR_KEEP_THE_HOST_ROUTE_NEGATIVE",
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
