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
        "R10": load_json("fundamental_action_reconstruction/generated/r10_legacy_control_to_current_pair_chart_reduction_packet_for_kobs_summary.json"),
        "P10": load_json("fundamental_action_reconstruction/generated/p10_existing_kernel_feedback_to_kobs_rerun_after_explicit_current_pair_chain_summary.json"),
        "P14": load_json("fundamental_action_reconstruction/generated/p14_existing_kernel_feedback_factorization_rerun_after_current_pair_chart_reduction_packet_summary.json"),
        "C10": load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json"),
        "H16": load_json("fundamental_action_reconstruction/generated/h16_operator_origin_partial_witness_audit_summary.json"),
        "H33": load_json("fundamental_action_reconstruction/generated/h33_pair1_selector_target_justification_audit.json"),
        "H34": load_json("fundamental_action_reconstruction/generated/h34_basis_covariance_target_independence_audit.json"),
    }

    checks_spec = [
        {
            "id": "r10_current_pair_chart_reduction_present",
            "actual": sources["R10"]["result"],
            "expected": "explicit_current_pair_chart_reduction_present_but_chart_scoped_only_and_not_a_strict_selector_target_justification",
            "meaning": "the chosen explicit current-pair chart reduction is present",
        },
        {
            "id": "p10_full_block_still_present",
            "actual": sources["P10"]["status"],
            "expected": "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "meaning": "the explicit current-pair H3 block is still present on the route",
        },
        {
            "id": "p14_factorization_route_still_negative",
            "actual": sources["P14"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_TO_EXPLICIT_CHAIN_FACTORIZATION_ROUTE_AFTER_CURRENT_PAIR_CHART_REDUCTION_PACKET",
            "meaning": "the factorization route remains noncomputable after current-pair chart reduction",
        },
        {
            "id": "c10_host_to_concrete_block_identification_not_shown",
            "actual": sources["C10"]["result"]["host_to_concrete_psi_block_identification"],
            "expected": "not_shown",
            "meaning": "no coefficient-level legacy-host to concrete-Psi-block identification is exported",
        },
        {
            "id": "h16_only_partial_witnesses_present",
            "actual": sources["H16"]["result"],
            "expected": "admissible_values_have_asymmetric_partial_witnesses_only",
            "meaning": "only partial operator-origin witnesses remain on the legacy side",
        },
        {
            "id": "h33_pair1_local_chart_only",
            "actual": sources["H33"]["result"],
            "expected": "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "meaning": "pair1 remains only a local chart",
        },
        {
            "id": "h34_no_covariance_argument",
            "actual": sources["H34"]["status"],
            "expected": "PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT",
            "meaning": "no basis-covariant or target-independent promotion is present",
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
        / "n17_current_existing_kernel_feedback_factorization_obstruction_after_current_pair_chart_reduction_packet_theorem_summary.json"
    )

    if mismatches:
        summary = {
            "step": "N17",
            "status": "N17_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FACTORIZATION_ROUTE_FRONTIER_AFTER_CURRENT_PAIR_CHART_REDUCTION_PACKET",
            "goal": "Check whether the current repo identifies existing kernel feedback with the explicit selector-facing H3 chain after current-pair chart reduction.",
            "scope": "current_existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_R10_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected factorization-route frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_FACTORIZATION_ROUTE_FRONTIER_AFTER_R10_BEFORE_CLAIMING_N17",
        }
    else:
        summary = {
            "step": "N17",
            "status": "N17_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_FACTORIZATION_OBSTRUCTION_AFTER_CURRENT_PAIR_CHART_REDUCTION_PACKET_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after materializing the chosen current-pair chart reduction, the current repo still does not identify existing kernel feedback with the explicit selector-facing H3 chain.",
            "scope": "current_existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_R10_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "chosen_current_pair_chart_reduction_present": True,
                "full_current_pair_h3_block_present": True,
                "existing_kernel_feedback_identified_with_explicit_chain": False,
            },
            "missing_structure_classes": [
                "intertwiner_or_equality_witness_identifying_the_chart_reduced_legacy_object_with_the_computed_current_pair_H3_block"
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future factorization routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "IMPLEMENT_THE_SINGLE_REMAINING_INTERTWINER_OR_EQUALITY_WITNESS_AND_RERUN_P14",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
