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
        "P15": load_json("fundamental_action_reconstruction/generated/p15_existing_kernel_feedback_intertwiner_equality_witness_probe_summary.json"),
        "C10": load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json"),
        "C15": load_json("fundamental_action_reconstruction/generated/c15_control_only_pullback_submatrix_packet_summary.json"),
        "H15": load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json"),
        "H18": load_json("fundamental_action_reconstruction/generated/h18_composite_route_a_provenance_binding_instance_summary.json"),
        "O2": load_json("fundamental_action_reconstruction/generated/o2_exported_composite_a1_ext_instance.json"),
    }

    checks_spec = [
        {
            "id": "r10_current_pair_chart_reduction_present",
            "actual": sources["R10"]["result"],
            "expected": "explicit_current_pair_chart_reduction_present_but_chart_scoped_only_and_not_a_strict_selector_target_justification",
            "meaning": "the chosen current-pair chart reduction is present",
        },
        {
            "id": "p10_current_pair_block_present",
            "actual": sources["P10"]["status"],
            "expected": "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "meaning": "the computed current-pair H3 block is present",
        },
        {
            "id": "p14_last_nominal_witness_blocker_present",
            "actual": sources["P14"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_TO_EXPLICIT_CHAIN_FACTORIZATION_ROUTE_AFTER_CURRENT_PAIR_CHART_REDUCTION_PACKET",
            "meaning": "P14 leaves the route at the single nominal witness blocker",
        },
        {
            "id": "p15_intertwiner_equality_route_negative",
            "actual": sources["P15"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_INTERTWINER_EQUALITY_WITNESS_ROUTE",
            "meaning": "the dedicated intertwiner/equality witness route remains noncomputable",
        },
        {
            "id": "c10_host_to_concrete_block_identification_not_shown",
            "actual": sources["C10"]["result"]["host_to_concrete_psi_block_identification"],
            "expected": "not_shown",
            "meaning": "the legacy host operator is still not identified with a concrete Psi-sector block",
        },
        {
            "id": "c15_coefficient_filled_control_pullback_present",
            "actual": sources["C15"]["result"]["coefficient_filled_M_control_present"],
            "expected": "yes",
            "meaning": "a coefficient-filled declared control pullback M_control is exported (declared scope; see P476 pointer in C15)",
        },
        {
            "id": "h15_existing_feedback_not_identified_with_kobs",
            "actual": sources["H15"]["result"],
            "expected": "existing_feedback_not_identified_with_kobs",
            "meaning": "existing kernel feedback still has no selector-facing export identified as its own reduction",
        },
        {
            "id": "h18_extension_lane_route_a_still_unevaluated",
            "actual": sources["H18"]["result"],
            "expected": "route_a_binding_completed_but_coefficients_unevaluated",
            "meaning": "the strongest extension-lane composite witness remains unevaluated",
        },
        {
            "id": "o2_composite_instance_still_unresolved",
            "actual": bool(sources["O2"]["computable_now"]),
            "expected": False,
            "meaning": "the exported composite A_1 instance is still not computable now",
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
        / "n18_current_existing_kernel_feedback_intertwiner_equality_witness_nonderivation_theorem_summary.json"
    )

    if mismatches:
        summary = {
            "step": "N18",
            "status": "N18_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_INTERTWINER_EQUALITY_WITNESS_FRONTIER",
            "goal": "Check whether the current repo exports an intertwiner/equality witness identifying the chart-reduced legacy object with the computed current-pair H3 block.",
            "scope": "current_existing_kernel_feedback_intertwiner_equality_witness_route_after_R10_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected witness-route frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_INTERTWINER_EQUALITY_WITNESS_FRONTIER_BEFORE_CLAIMING_N18",
        }
    else:
        summary = {
            "step": "N18",
            "status": "N18_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_INTERTWINER_EQUALITY_WITNESS_NONDERIVATION_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after reaching the chosen current-pair chart, the current repo still does not export an intertwiner/equality witness identifying existing kernel feedback with the computed current-pair H3 block.",
            "scope": "current_existing_kernel_feedback_intertwiner_equality_witness_route_after_R10_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "chosen_current_pair_chart_reduction_present": True,
                "computed_current_pair_h3_block_present": True,
                "legacy_chart_reduced_operator_object_present": False,
                "intertwiner_or_equality_witness_present": False,
            },
            "missing_structure_classes": [
                "explicit_coefficient_filled_legacy_chart_reduced_operator_object_on_the_chosen_current_pair_chart_pair1_or_equivalent_actual_target",
                "intertwiner_or_equality_witness_identifying_that_legacy_chart_reduced_operator_object_with_the_computed_current_pair_H3_block",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future factorization routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "EXPORT_THE_LEGACY_CHART_REDUCED_OPERATOR_OBJECT_FIRST_OR_KEEP_THE_EXISTING_KERNEL_FEEDBACK_FACTORIZATION_ROUTE_NEGATIVE",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
