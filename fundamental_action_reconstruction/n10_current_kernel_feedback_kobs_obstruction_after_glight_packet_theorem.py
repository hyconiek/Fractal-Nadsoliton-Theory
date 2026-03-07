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
        "R3": load_json("fundamental_action_reconstruction/generated/r3_minimal_internal_light_propagator_packet_for_kobs_summary.json"),
        "R3_artifact": load_json("fundamental_action_reconstruction/generated/r3_minimal_internal_light_propagator_packet_for_kobs.json"),
        "H14": load_json("fundamental_action_reconstruction/generated/h14_existing_kernel_feedback_vs_new_kobs_separation_audit_summary.json"),
        "H15": load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json"),
        "H29": load_json("fundamental_action_reconstruction/generated/h29_wave_retardation_proxy_selector_reduction_audit_summary.json"),
        "P7": load_json("fundamental_action_reconstruction/generated/p7_existing_kernel_feedback_to_kobs_rerun_after_glight_packet_summary.json"),
    }

    checks_spec = [
        {
            "id": "p7_route_negative",
            "actual": sources["P7"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_KERNEL_FEEDBACK_TO_KOBS_ROUTE_AFTER_GLIGHT_PACKET",
            "meaning": "the rerun probe still does not instantiate a selector-facing K_obs",
        },
        {
            "id": "r3_explicit_glight_packet_present",
            "actual": sources["R3"]["result"],
            "expected": "explicit_internal_light_propagator_packet_present_but_eigenchannel_only_and_not_yet_factorized_from_existing_kernel_feedback",
            "meaning": "an explicit finite G_light packet exists",
        },
        {
            "id": "r3_factorization_still_absent",
            "actual": sources["R3_artifact"]["factorization_status"],
            "expected": "not_identified_with_existing_kernel_feedback",
            "meaning": "the new packet is not yet a factorization of current kernel feedback",
        },
        {
            "id": "r3_matrix_no_orientation_anchor",
            "actual": sources["R3_artifact"]["uses_orientation_anchor_in_matrix"],
            "expected": False,
            "meaning": "the explicit G_light matrix does not use psi0",
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

    out = ROOT / "generated" / "n10_current_kernel_feedback_kobs_obstruction_after_glight_packet_theorem_summary.json"

    if mismatches:
        summary = {
            "step": "N10",
            "status": "N10_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_UPDATED_KOBS_ROUTE_FRONTIER",
            "goal": "Check whether the updated kernel-feedback-to-K_obs route instantiates a selector-facing K_obs after adding an explicit G_light packet.",
            "scope": "current_kernel_feedback_to_kobs_route_after_R3_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected updated K_obs frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_UPDATED_KOBS_ROUTE_FRONTIER_BEFORE_CLAIMING_N10",
        }
    else:
        summary = {
            "step": "N10",
            "status": "N10_DISCHARGED_CURRENT_KERNEL_FEEDBACK_KOBS_OBSTRUCTION_AFTER_GLIGHT_PACKET_NO_FALSE_PASS",
            "goal": "Discharge an updated route-specific theorem: even after adding an explicit G_light packet, the current kernel-feedback route does not yet instantiate a selector-facing K_obs.",
            "scope": "current_kernel_feedback_to_kobs_route_after_R3_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "explicit_internal_light_propagator_packet_present": True,
                "selector_facing_kobs_instantiated": False,
                "explicit_glight_packet_is_not_factorization_discharge": True,
                "preoriented_proxies_are_not_internal_anchor_generation": True,
            },
            "missing_structure_classes": [
                "explicit_emission_map_E_from_M_pair_to_L_int",
                "explicit_light_to_matter_response_map_R_mat_from_L_int_to_Q_mat",
                "explicit_observer_readout_operator_O_obs_on_Q_mat",
                "equivalence_or_factorization_map_from_existing_kernel_feedback_and_R2_parameter_packet_to_H3_operator_chain",
                "selector_sector_projected_2x2_block_export_on_an_actual_pair",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future K_obs factorizations are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "ADD_ONE_REMAINING_OPERATOR_CHAIN_OBJECT_OR_THE_FACTORIZATION_MAP_AND_RERUN_P7",
        }

    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
