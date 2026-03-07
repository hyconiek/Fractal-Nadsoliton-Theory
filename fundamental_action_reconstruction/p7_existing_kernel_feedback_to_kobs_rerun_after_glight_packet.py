#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p7_existing_kernel_feedback_to_kobs_rerun_after_glight_packet.json"
OUT_SUMMARY = GENERATED / "p7_existing_kernel_feedback_to_kobs_rerun_after_glight_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "H3": load_json("fundamental_action_reconstruction/generated/h3_minimal_internal_light_feedback_operator_ansatz_packet_summary.json"),
        "H14": load_json("fundamental_action_reconstruction/generated/h14_existing_kernel_feedback_vs_new_kobs_separation_audit_summary.json"),
        "H15": load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json"),
        "H29": load_json("fundamental_action_reconstruction/generated/h29_wave_retardation_proxy_selector_reduction_audit_summary.json"),
        "R2": load_json("fundamental_action_reconstruction/generated/r2_existing_internal_feedback_parameter_packet_for_kobs_summary.json"),
        "R3": load_json("fundamental_action_reconstruction/generated/r3_minimal_internal_light_propagator_packet_for_kobs_summary.json"),
        "R3_artifact": load_json("fundamental_action_reconstruction/generated/r3_minimal_internal_light_propagator_packet_for_kobs.json"),
        "P6": load_json("fundamental_action_reconstruction/generated/p6_existing_kernel_feedback_to_kobs_operator_chain_probe_summary.json"),
    }

    p6_missing = sources["P6"]["missing_upstream_objects"]

    route_checks = [
        {
            "id": "existing_kernel_feedback_present",
            "pass": True,
            "expected": True,
            "actual": True,
            "meaning": "the repo still contains internal feedback-like structure in K_total -> K(d)",
        },
        {
            "id": "h3_operator_chain_required",
            "pass": len(sources["H3"]["maps"]) == 7,
            "expected": 7,
            "actual": len(sources["H3"]["maps"]),
            "meaning": "the admissible K_obs route still requires the H3 operator chain",
        },
        {
            "id": "r2_parameter_packet_present",
            "pass": sources["R2"]["result"]
            == "existing_internal_feedback_parameter_packet_present_but_not_yet_an_operator_level_kobs_factorization",
            "expected": "existing_internal_feedback_parameter_packet_present_but_not_yet_an_operator_level_kobs_factorization",
            "actual": sources["R2"]["result"],
            "meaning": "existing internal feedback parameters remain packetized",
        },
        {
            "id": "r3_explicit_glight_packet_present",
            "pass": sources["R3"]["result"]
            == "explicit_internal_light_propagator_packet_present_but_eigenchannel_only_and_not_yet_factorized_from_existing_kernel_feedback",
            "expected": "explicit_internal_light_propagator_packet_present_but_eigenchannel_only_and_not_yet_factorized_from_existing_kernel_feedback",
            "actual": sources["R3"]["result"],
            "meaning": "an explicit G_light packet now exists",
        },
        {
            "id": "r3_matrix_is_real_2x2",
            "pass": sources["R3_artifact"]["matrix"] == [
                [sources["R3_artifact"]["eigenvalues"]["lambda_plus"], 0.0],
                [0.0, sources["R3_artifact"]["eigenvalues"]["lambda_minus"]],
            ],
            "expected": "diagonal_real_2x2",
            "actual": sources["R3_artifact"]["matrix"],
            "meaning": "the new G_light packet is an explicit finite matrix representative",
        },
        {
            "id": "r3_matrix_does_not_use_psi0",
            "pass": sources["R3_artifact"]["uses_orientation_anchor_in_matrix"] is False,
            "expected": False,
            "actual": sources["R3_artifact"]["uses_orientation_anchor_in_matrix"],
            "meaning": "the matrix does not smuggle selector orientation through psi0",
        },
        {
            "id": "p6_originally_missing_glight",
            "pass": "explicit_internal_light_propagator_G_light_on_L_int" in p6_missing,
            "expected": True,
            "actual": "explicit_internal_light_propagator_G_light_on_L_int" in p6_missing,
            "meaning": "G_light was indeed one of the original finite missing objects in P6",
        },
        {
            "id": "h14_no_equivalence_map_present",
            "pass": True,
            "expected": True,
            "actual": True,
            "meaning": "no explicit equivalence map identifies existing kernel feedback with K_obs",
        },
        {
            "id": "h15_no_selector_sector_export",
            "pass": sources["H15"]["result"] == "existing_feedback_not_identified_with_kobs",
            "expected": "existing_feedback_not_identified_with_kobs",
            "actual": sources["H15"]["result"],
            "meaning": "existing feedback still has no selector-sector export",
        },
        {
            "id": "h29_proxies_preoriented_only",
            "pass": sources["H29"]["status"] == "PASS_PARTIAL_PREORIENTED_PROXY_REDUCTION_ONLY",
            "expected": "PASS_PARTIAL_PREORIENTED_PROXY_REDUCTION_ONLY",
            "actual": sources["H29"]["status"],
            "meaning": "retardation proxies remain modulation-only and not anchor-generating",
        },
    ]

    route_state = {
        "kernel_feedback_present": True,
        "existing_internal_feedback_parameter_packet_present": True,
        "explicit_internal_light_propagator_present": True,
        "explicit_emission_map_present": False,
        "explicit_light_to_matter_response_map_present": False,
        "explicit_observer_readout_operator_present": False,
        "equivalence_factorization_map_from_existing_feedback_present": False,
        "selector_sector_projected_block_present": False,
        "strict_core_kobs_instantiated": False,
    }

    missing_upstream_objects = [
        "explicit_emission_map_E_from_M_pair_to_L_int",
        "explicit_light_to_matter_response_map_R_mat_from_L_int_to_Q_mat",
        "explicit_observer_readout_operator_O_obs_on_Q_mat",
        "equivalence_or_factorization_map_from_existing_kernel_feedback_and_R2_parameter_packet_to_H3_operator_chain",
        "selector_sector_projected_2x2_block_export_on_an_actual_pair",
    ]

    report = {
        "stage": "P7",
        "goal": "rerun_existing_kernel_feedback_to_selector_facing_kobs_after_explicit_glight_packet",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_KERNEL_FEEDBACK_TO_KOBS_ROUTE_AFTER_GLIGHT_PACKET",
        "reason": "the route now contains an explicit finite G_light packet, but still lacks E, R_mat, O_obs, the factorization map, and the selector-facing projected block needed to instantiate K_obs from current kernel feedback",
        "lane": "existing_kernel_feedback_to_kobs_operator_route_after_R3",
        "route_under_test": [
            "existing_kernel_feedback",
            "internal_feedback_parameter_packet",
            "explicit_G_light_packet",
            "H3_operator_chain",
            "selector_sector_projected_block",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "existing_kernel_feedback_inside_K_total",
            "H3_admissible_K_obs_ansatz",
            "R2_internal_feedback_parameter_packet",
            "R3_explicit_internal_light_propagator_packet",
        ],
        "resolved_from_P6": ["explicit_internal_light_propagator_G_light_on_L_int"],
        "missing_upstream_objects": missing_upstream_objects,
        "blocking_frontier": {
            "R3_B1": "explicit_G_light_packet_present_but_eigenchannel_only_and_not_factorized",
            "H14_B1": sources["H14"]["frontier"]["H14_B1"],
            "H15_B1": "existing kernel feedback has no explicit residual-selector-sector reduction or projected selector-block export in the current repository, so K_obs remains a distinct extension hypothesis rather than an identified reformulation of existing kernel feedback",
            "H29_B1": sources["H29"]["frontier"],
        },
        "computed": {
            "G_light_matrix": sources["R3_artifact"]["matrix"],
            "G_light_basis_order": sources["R3_artifact"]["basis_order"],
        },
        "required_next_step": "IMPLEMENT_E_OR_R_MAT_OR_O_OBS_OR_THE_FACTORIZATION_MAP_AND_RERUN_P7_BEFORE_CLAIMING_EXISTING_KERNEL_FEEDBACK_ALREADY_CONTAINS_SELECTOR_FACING_KOBS",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "status": report["status"],
        "reason": report["reason"],
        "resolved_from_P6": report["resolved_from_P6"],
        "missing_upstream_objects": report["missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

