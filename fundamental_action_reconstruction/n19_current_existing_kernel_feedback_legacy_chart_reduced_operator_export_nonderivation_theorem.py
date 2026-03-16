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
        "R8": load_json("fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs_summary.json"),
        "R9": load_json("fundamental_action_reconstruction/generated/r9_existing_kernel_feedback_host_to_control_pushforward_packet_for_kobs_summary.json"),
        "R10": load_json("fundamental_action_reconstruction/generated/r10_legacy_control_to_current_pair_chart_reduction_packet_for_kobs_summary.json"),
        "P15": load_json("fundamental_action_reconstruction/generated/p15_existing_kernel_feedback_intertwiner_equality_witness_probe_summary.json"),
        "P16": load_json("fundamental_action_reconstruction/generated/p16_existing_kernel_feedback_legacy_chart_reduced_operator_export_probe_summary.json"),
        "C10": load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json"),
        "C11": load_json("fundamental_action_reconstruction/generated/c11_psi_sector_block_extraction_audit_summary.json"),
        "C15": load_json("fundamental_action_reconstruction/generated/c15_control_only_pullback_submatrix_packet_summary.json"),
        "C20": load_json("fundamental_action_reconstruction/generated/c20_finite_materialization_recipe_audit_summary.json"),
    }

    checks_spec = [
        {
            "id": "r8_host_operator_carrier_present",
            "actual": sources["R8"]["result"],
            "expected": "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
            "meaning": "the legacy host carrier is present",
        },
        {
            "id": "r9_typed_pushforward_present",
            "actual": sources["R9"]["result"],
            "expected": "typed_host_to_control_pushforward_present_but_control_scope_only_and_not_yet_reduced_to_the_explicit_current_pair_chain",
            "meaning": "the typed host-to-control pushforward is present",
        },
        {
            "id": "r10_current_pair_chart_reduction_present",
            "actual": sources["R10"]["result"],
            "expected": "explicit_current_pair_chart_reduction_present_but_chart_scoped_only_and_not_a_strict_selector_target_justification",
            "meaning": "the chosen current-pair chart reduction is present",
        },
        {
            "id": "p15_legacy_chart_reduced_operator_object_was_the_missing_blocker",
            "actual": sources["P15"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_INTERTWINER_EQUALITY_WITNESS_ROUTE",
            "meaning": "P15 localized the missing legacy chart-reduced operator object upstream of the witness",
        },
        {
            "id": "p16_legacy_operator_export_route_negative",
            "actual": sources["P16"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_LEGACY_CHART_REDUCED_OPERATOR_EXPORT_ROUTE",
            "meaning": "the dedicated legacy chart-reduced operator export route remains noncomputable",
        },
        {
            "id": "c10_host_to_concrete_psi_block_identification_not_shown",
            "actual": sources["C10"]["result"]["host_to_concrete_psi_block_identification"],
            "expected": "not_shown",
            "meaning": "the host is still not identified with a concrete Psi-sector block",
        },
        {
            "id": "c11_concrete_psi_block_not_extracted",
            "actual": sources["C11"]["result"]["concrete_block_extracted"],
            "expected": "not_shown",
            "meaning": "no concrete Psi-sector block has been extracted",
        },
        {
            "id": "c15_coefficient_filled_control_pullback_present",
            "actual": sources["C15"]["result"]["coefficient_filled_M_control_present"],
            "expected": "yes",
            "meaning": "a coefficient-filled declared control pullback M_control is exported (declared scope; see P476 pointer in C15)",
        },
        {
            "id": "c20_executed_serialization_run_not_shown",
            "actual": sources["C20"]["persisted_outputs"]["persisted_12_row_serialization_run_present"],
            "expected": False,
            "meaning": "the finite materialization recipe still has no executed persisted 12-row serialization run",
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
        / "n19_current_existing_kernel_feedback_legacy_chart_reduced_operator_export_nonderivation_theorem_summary.json"
    )

    if mismatches:
        summary = {
            "step": "N19",
            "status": "N19_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_OPERATOR_EXPORT_FRONTIER",
            "goal": "Check whether the current repo exports a coefficient-filled legacy chart-reduced operator object on pair1 from existing kernel feedback.",
            "scope": "current_existing_kernel_feedback_legacy_chart_reduced_operator_export_route_after_P15_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected legacy-operator-export frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_LEGACY_OPERATOR_EXPORT_FRONTIER_BEFORE_CLAIMING_N19",
        }
    else:
        summary = {
            "step": "N19",
            "status": "N19_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_LEGACY_CHART_REDUCED_OPERATOR_EXPORT_NONDERIVATION_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after host-carrier packetization, typed pushforward, and current-pair chart reduction, the current repo still does not export a coefficient-filled legacy chart-reduced operator object on pair1 from existing kernel feedback.",
            "scope": "current_existing_kernel_feedback_legacy_chart_reduced_operator_export_route_after_P15_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "host_operator_carrier_present": True,
                "typed_host_to_control_pushforward_present": True,
                "chosen_current_pair_chart_reduction_present": True,
                "coefficient_filled_legacy_chart_reduced_operator_present": False,
            },
            "missing_structure_classes": [
                "host_to_concrete_Psi_sector_quadratic_block_identification_for_the_existing_kernel_feedback_host_operator"
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future factorization routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "EITHER_EXPORT_A_HOST_TO_CANONICAL_PSI_BLOCK_MATCHING_WITNESS_IDENTIFYING_QW2186_WITH_THE_CANONICAL_CARRIER_OR_KEEP_THE_FACTORING_ROUTE_NEGATIVE",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
