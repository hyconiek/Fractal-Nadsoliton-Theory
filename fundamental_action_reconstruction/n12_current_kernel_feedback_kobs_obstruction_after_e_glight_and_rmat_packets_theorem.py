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
        "R4": load_json("fundamental_action_reconstruction/generated/r4_local_chart_emission_map_packet_for_kobs_summary.json"),
        "R4_artifact": load_json("fundamental_action_reconstruction/generated/r4_local_chart_emission_map_packet_for_kobs.json"),
        "R5": load_json("fundamental_action_reconstruction/generated/r5_minimal_light_to_matter_response_packet_for_kobs_summary.json"),
        "R5_artifact": load_json("fundamental_action_reconstruction/generated/r5_minimal_light_to_matter_response_packet_for_kobs.json"),
        "H15": load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json"),
        "H29": load_json("fundamental_action_reconstruction/generated/h29_wave_retardation_proxy_selector_reduction_audit_summary.json"),
        "H33": load_json("fundamental_action_reconstruction/generated/h33_pair1_selector_target_justification_audit.json"),
        "H34": load_json("fundamental_action_reconstruction/generated/h34_basis_covariance_target_independence_audit.json"),
        "H35": load_json("fundamental_action_reconstruction/generated/h35_pair1_axis_selection_audit.json"),
        "H36": load_json("fundamental_action_reconstruction/generated/h36_directed_axis_orientation_audit_summary.json"),
        "H37": load_json("fundamental_action_reconstruction/generated/h37_sign_distinction_state_audit_summary.json"),
        "P9": load_json("fundamental_action_reconstruction/generated/p9_existing_kernel_feedback_to_kobs_rerun_after_e_glight_and_rmat_packets_summary.json"),
    }

    checks_spec = [
        {
            "id": "p9_route_negative",
            "actual": sources["P9"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_KERNEL_FEEDBACK_TO_KOBS_ROUTE_AFTER_E_GLIGHT_AND_RMAT_PACKETS",
            "meaning": "the rerun probe still does not instantiate a selector-facing K_obs",
        },
        {
            "id": "r4_explicit_e_packet_present",
            "actual": sources["R4"]["result"],
            "expected": "explicit_local_chart_emission_map_packet_present_but_preoriented_and_not_yet_factorized_from_existing_kernel_feedback",
            "meaning": "an explicit local-chart emission packet exists",
        },
        {
            "id": "r4_factorization_still_absent",
            "actual": sources["R4_artifact"]["factorization_status"],
            "expected": "not_identified_with_existing_kernel_feedback",
            "meaning": "the E packet is still not a factorization of current kernel feedback",
        },
        {
            "id": "r5_explicit_rmat_packet_present",
            "actual": sources["R5"]["result"],
            "expected": "explicit_light_to_matter_response_packet_present_but_current_pair_scoped_and_not_yet_factorized_from_existing_kernel_feedback",
            "meaning": "an explicit current-pair light-to-matter response packet exists",
        },
        {
            "id": "r5_factorization_still_absent",
            "actual": sources["R5_artifact"]["factorization_status"],
            "expected": "not_identified_with_existing_kernel_feedback",
            "meaning": "the new R_mat packet is still not a factorization of current kernel feedback",
        },
        {
            "id": "r5_pair_scope_only",
            "actual": bool(sources["R5_artifact"]["current_pair_scope_only"]),
            "expected": True,
            "meaning": "the new R_mat packet remains only current-pair scoped",
        },
        {
            "id": "h33_local_chart_only",
            "actual": sources["H33"]["status"],
            "expected": "PASS_PARTIAL_LOCAL_CHART_ONLY",
            "meaning": "pair1 remains only a local chart",
        },
        {
            "id": "h34_no_covariance_argument",
            "actual": sources["H34"]["status"],
            "expected": "PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT",
            "meaning": "basis-covariance remains undischarged",
        },
        {
            "id": "h35_no_axis_selection_argument",
            "actual": sources["H35"]["status"],
            "expected": "PASS_PARTIAL_NO_STRICT_AXIS_SELECTION_ARGUMENT",
            "meaning": "strict axis selection remains absent",
        },
        {
            "id": "h36_directed_axis_orientation_present_premise_based",
            "actual": sources["H36"]["status"],
            "expected": "PASS_DIRECTED_AXIS_ORIENTATION_PRESENT__PREMISE_BASED_T164",
            "meaning": "directed-axis orientation mechanism is now exported (premise-based), but this does not instantiate selector-facing K_obs from existing kernel feedback",
        },
        {
            "id": "h37_sign_distinction_state_exported_premise_based",
            "actual": sources["H37"]["status"],
            "expected": "PASS_STRICT_DIRECTED_SIGN_SENSITIVE_DISTINCTION_OBSERVABLE_EXPORTED_AND_GLOBAL_DIRECTED_SELECTOR_STATE_OBJECT_EXPORTED__PREMISE_BASED_T164",
            "meaning": "a sign-sensitive directed selector-state layer is now exported (premise-based), but this does not instantiate selector-facing K_obs from existing kernel feedback",
        },
        {
            "id": "h15_selector_sector_export_absent",
            "actual": sources["H15"]["result"],
            "expected": "existing_feedback_not_identified_with_kobs",
            "meaning": "selector-sector export remains absent",
        },
        {
            "id": "h29_proxies_preoriented_only",
            "actual": sources["H29"]["status"],
            "expected": "PASS_PARTIAL_PREORIENTED_PROXY_REDUCTION_ONLY",
            "meaning": "retardation proxies remain modulation-only and not anchor-generating",
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

    out = ROOT / "generated" / "n12_current_kernel_feedback_kobs_obstruction_after_e_glight_and_rmat_packets_theorem_summary.json"

    if mismatches:
        summary = {
            "step": "N12",
            "status": "N12_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_UPDATED_KOBS_ROUTE_FRONTIER",
            "goal": "Check whether the updated kernel-feedback-to-K_obs route instantiates a selector-facing K_obs after adding explicit E, G_light, and R_mat packets.",
            "scope": "current_kernel_feedback_to_kobs_route_after_R5_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected updated K_obs frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_UPDATED_KOBS_ROUTE_FRONTIER_BEFORE_CLAIMING_N12",
        }
    else:
        summary = {
            "step": "N12",
            "status": "N12_DISCHARGED_CURRENT_KERNEL_FEEDBACK_KOBS_OBSTRUCTION_AFTER_E_GLIGHT_AND_RMAT_PACKETS_NO_FALSE_PASS",
            "goal": "Discharge an updated route-specific theorem: even after adding explicit E, G_light, and R_mat packets, the current kernel-feedback route does not yet instantiate a selector-facing K_obs.",
            "scope": "current_kernel_feedback_to_kobs_route_after_R5_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "explicit_emission_map_packet_present": True,
                "explicit_internal_light_propagator_packet_present": True,
                "explicit_light_to_matter_response_packet_present": True,
                "emission_packet_remains_local_chart_preoriented": True,
                "matter_response_packet_remains_current_pair_scoped": True,
                "selector_facing_kobs_instantiated": False,
            },
            "missing_structure_classes": [
                "explicit_observer_readout_operator_O_obs_on_Q_mat",
                "equivalence_or_factorization_map_from_existing_kernel_feedback_and_R2_parameter_packet_to_H3_operator_chain",
                "full_H3_selector_sector_projected_2x2_block_export_on_an_actual_pair",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future K_obs factorizations are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "ADD_O_OBS_OR_THE_FACTORIZATION_MAP_AND_RERUN_P9",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
