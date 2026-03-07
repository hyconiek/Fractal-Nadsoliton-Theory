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
        "H14": load_json("fundamental_action_reconstruction/generated/h14_existing_kernel_feedback_vs_new_kobs_separation_audit_summary.json"),
        "H15": load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json"),
        "H29": load_json("fundamental_action_reconstruction/generated/h29_wave_retardation_proxy_selector_reduction_audit_summary.json"),
        "R2": load_json("fundamental_action_reconstruction/generated/r2_existing_internal_feedback_parameter_packet_for_kobs_summary.json"),
        "P6": load_json("fundamental_action_reconstruction/generated/p6_existing_kernel_feedback_to_kobs_operator_chain_probe_summary.json"),
    }

    checks_spec = [
        {
            "id": "p6_route_negative",
            "actual": sources["P6"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_KERNEL_FEEDBACK_TO_KOBS_ROUTE",
            "meaning": "the route still does not instantiate selector-facing K_obs",
        },
        {
            "id": "r2_parameter_packet_present",
            "actual": sources["R2"]["result"],
            "expected": "existing_internal_feedback_parameter_packet_present_but_not_yet_an_operator_level_kobs_factorization",
            "meaning": "internal feedback parameters are now collected in one packet",
        },
        {
            "id": "h14_non_identification_active",
            "actual": sources["H14"]["status"],
            "expected": "PASS_PARTIAL_EXISTING_FEEDBACK_RECOGNIZED_BUT_NOT_IDENTIFIED_WITH_KOBS",
            "meaning": "existing kernel feedback is still not identified with K_obs",
        },
        {
            "id": "h15_no_selector_sector_export",
            "actual": sources["H15"]["result"],
            "expected": "existing_feedback_not_identified_with_kobs",
            "meaning": "existing feedback is still not exported to the selector sector",
        },
        {
            "id": "h29_preoriented_proxy_only",
            "actual": sources["H29"]["status"],
            "expected": "PASS_PARTIAL_PREORIENTED_PROXY_REDUCTION_ONLY",
            "meaning": "old proxies remain preoriented only and do not generate an internal anchor",
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
            "step": "N9",
            "status": "N9_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_KERNEL_FEEDBACK_TO_KOBS_FRONTIER",
            "goal": "Check whether current kernel feedback already instantiates a selector-facing K_obs operator.",
            "scope": "current_kernel_feedback_to_kobs_route_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected kernel-feedback-to-K_obs frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_KERNEL_FEEDBACK_TO_KOBS_FRONTIER_BEFORE_CLAIMING_N9",
        }
    else:
        summary = {
            "step": "N9",
            "status": "N9_DISCHARGED_CURRENT_KERNEL_FEEDBACK_DOES_NOT_YET_INSTANTIATE_SELECTOR_FACING_KOBS_NO_FALSE_PASS",
            "goal": "Discharge a current-route theorem: existing kernel feedback does not yet instantiate a selector-facing K_obs operator.",
            "scope": "current_kernel_feedback_to_kobs_route_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "existing_feedback_parameters_present": True,
                "selector_facing_kobs_instantiated": False,
                "parameter_packet_is_not_operator_factorization": True,
                "preoriented_proxies_are_not_internal_anchor_generation": True,
            },
            "missing_structure_classes": [
                "explicit_emission_map_E",
                "explicit_internal_light_propagator_G_light",
                "explicit_light_to_matter_response_map_R_mat",
                "explicit_observer_readout_operator_O_obs",
                "equivalence_or_factorization_map_from_existing_feedback_to_H3_chain",
                "selector_sector_projected_block_export",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future kernel-feedback factorization into K_obs is impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "ADD_ONE_OPERATOR_CHAIN_OBJECT_OR_FACTORIZATION_MAP_AND_RERUN_P6_BEFORE_PROMOTING_KOBS_BEYOND_LIVE_EXTENSION_HYPOTHESIS",
        }

    out = ROOT / "generated" / "n9_current_kernel_feedback_does_not_yet_instantiate_selector_facing_kobs_theorem_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
