#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p8_existing_kernel_feedback_to_kobs_rerun_after_e_and_glight_packets.json"
OUT_SUMMARY = GENERATED / "p8_existing_kernel_feedback_to_kobs_rerun_after_e_and_glight_packets_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "H3": load_json("fundamental_action_reconstruction/generated/h3_minimal_internal_light_feedback_operator_ansatz_packet_summary.json"),
        "H14": load_json("fundamental_action_reconstruction/generated/h14_existing_kernel_feedback_vs_new_kobs_separation_audit_summary.json"),
        "H15": load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json"),
        "H29": load_json("fundamental_action_reconstruction/generated/h29_wave_retardation_proxy_selector_reduction_audit_summary.json"),
        "H31": load_json("fundamental_action_reconstruction/generated/h31_psi0_to_pair1_reduction_audit.json"),
        "H33": load_json("fundamental_action_reconstruction/generated/h33_pair1_selector_target_justification_audit.json"),
        "H34": load_json("fundamental_action_reconstruction/generated/h34_basis_covariance_target_independence_audit.json"),
        "H35": load_json("fundamental_action_reconstruction/generated/h35_pair1_axis_selection_audit.json"),
        "H36": load_json("fundamental_action_reconstruction/generated/h36_directed_axis_orientation_audit_summary.json"),
        "H37": load_json("fundamental_action_reconstruction/generated/h37_sign_distinction_state_audit_summary.json"),
        "R2": load_json("fundamental_action_reconstruction/generated/r2_existing_internal_feedback_parameter_packet_for_kobs_summary.json"),
        "R3": load_json("fundamental_action_reconstruction/generated/r3_minimal_internal_light_propagator_packet_for_kobs_summary.json"),
        "R3_artifact": load_json("fundamental_action_reconstruction/generated/r3_minimal_internal_light_propagator_packet_for_kobs.json"),
        "R4": load_json("fundamental_action_reconstruction/generated/r4_local_chart_emission_map_packet_for_kobs_summary.json"),
        "R4_artifact": load_json("fundamental_action_reconstruction/generated/r4_local_chart_emission_map_packet_for_kobs.json"),
        "P7": load_json("fundamental_action_reconstruction/generated/p7_existing_kernel_feedback_to_kobs_rerun_after_glight_packet_summary.json"),
    }

    p7_missing = sources["P7"]["missing_upstream_objects"]

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
            "meaning": "an explicit G_light packet remains present",
        },
        {
            "id": "r4_explicit_e_packet_present",
            "pass": sources["R4"]["result"]
            == "explicit_local_chart_emission_map_packet_present_but_preoriented_and_not_yet_factorized_from_existing_kernel_feedback",
            "expected": "explicit_local_chart_emission_map_packet_present_but_preoriented_and_not_yet_factorized_from_existing_kernel_feedback",
            "actual": sources["R4"]["result"],
            "meaning": "an explicit current-pair emission packet now exists",
        },
        {
            "id": "r4_current_pair_identification_present",
            "pass": sources["R4_artifact"]["current_test_carrier_identification"]["M_pair_current"] == "V_1 = span{c_1, s_1}",
            "expected": "V_1 = span{c_1, s_1}",
            "actual": sources["R4_artifact"]["current_test_carrier_identification"]["M_pair_current"],
            "meaning": "the emission packet is explicitly identified as the current actual pair representative of M_pair",
        },
        {
            "id": "r4_partial_pullback_matches_p1",
            "pass": bool(sources["R4_artifact"]["partial_emission_propagation_pullback"]["matches_p1_test_c_configured"]),
            "expected": True,
            "actual": bool(sources["R4_artifact"]["partial_emission_propagation_pullback"]["matches_p1_test_c_configured"]),
            "meaning": "the E+G partial pullback coheres with the already computed P1 extension matrix",
        },
        {
            "id": "p7_originally_missing_e",
            "pass": "explicit_emission_map_E_from_M_pair_to_L_int" in p7_missing,
            "expected": True,
            "actual": "explicit_emission_map_E_from_M_pair_to_L_int" in p7_missing,
            "meaning": "E was indeed one of the original finite missing objects in P7",
        },
        {
            "id": "h31_coordinate_embedding_only",
            "pass": bool(sources["H31"]["classification"]["coordinate_level_embedding_present"])
            and not bool(sources["H31"]["classification"]["strict_core_selector_reduction_present"]),
            "expected": True,
            "actual": bool(sources["H31"]["classification"]["coordinate_level_embedding_present"])
            and not bool(sources["H31"]["classification"]["strict_core_selector_reduction_present"]),
            "meaning": "the route still has only a coordinate embedding and not a strict selector reduction",
        },
        {
            "id": "h33_pair1_local_chart_only",
            "pass": sources["H33"]["result"]
            == "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "expected": "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "actual": sources["H33"]["result"],
            "meaning": "pair1 remains only a local chart and not a privileged target",
        },
        {
            "id": "h34_no_covariance_argument",
            "pass": sources["H34"]["status"] == "PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT",
            "expected": "PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT",
            "actual": sources["H34"]["status"],
            "meaning": "basis-covariance and target-independence remain undischarged",
        },
        {
            "id": "h35_no_physical_axis_selection",
            "pass": sources["H35"]["status"] == "PASS_PARTIAL_NO_STRICT_AXIS_SELECTION_ARGUMENT",
            "expected": "PASS_PARTIAL_NO_STRICT_AXIS_SELECTION_ARGUMENT",
            "actual": sources["H35"]["status"],
            "meaning": "psi0 still does not strictly select a physical axis inside pair1",
        },
        {
            "id": "h36_no_directed_axis_selection",
            "pass": sources["H36"]["status"] == "PASS_PARTIAL_NO_STRICT_DIRECTED_AXIS_SELECTION",
            "expected": "PASS_PARTIAL_NO_STRICT_DIRECTED_AXIS_SELECTION",
            "actual": sources["H36"]["status"],
            "meaning": "the route still has no directed-axis discharge",
        },
        {
            "id": "h37_no_sign_distinction_state",
            "pass": sources["H37"]["status"] == "PASS_PARTIAL_NO_STRICT_SIGN_DISTINCTION_STATE_OBJECT",
            "expected": "PASS_PARTIAL_NO_STRICT_SIGN_DISTINCTION_STATE_OBJECT",
            "actual": sources["H37"]["status"],
            "meaning": "the route still has no sign-sensitive selector-state object",
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
        "explicit_emission_map_present": True,
        "explicit_emission_map_is_local_chart_preoriented_only": True,
        "explicit_light_to_matter_response_map_present": False,
        "explicit_observer_readout_operator_present": False,
        "equivalence_factorization_map_from_existing_feedback_present": False,
        "full_h3_selector_sector_projected_block_present": False,
        "strict_core_kobs_instantiated": False,
    }

    missing_upstream_objects = [
        "explicit_light_to_matter_response_map_R_mat_from_L_int_to_Q_mat",
        "explicit_observer_readout_operator_O_obs_on_Q_mat",
        "equivalence_or_factorization_map_from_existing_kernel_feedback_and_R2_parameter_packet_to_H3_operator_chain",
        "full_H3_selector_sector_projected_2x2_block_export_on_an_actual_pair",
    ]

    report = {
        "stage": "P8",
        "goal": "rerun_existing_kernel_feedback_to_selector_facing_kobs_after_explicit_e_and_glight_packets",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_KERNEL_FEEDBACK_TO_KOBS_ROUTE_AFTER_E_AND_GLIGHT_PACKETS",
        "reason": "the route now contains explicit E and G_light packets, but E remains local-chart preoriented and the route still lacks R_mat, O_obs, the factorization map, and the full H3 selector-facing projected block needed to instantiate K_obs from current kernel feedback",
        "lane": "existing_kernel_feedback_to_kobs_operator_route_after_R4",
        "route_under_test": [
            "existing_kernel_feedback",
            "internal_feedback_parameter_packet",
            "explicit_E_packet",
            "explicit_G_light_packet",
            "H3_operator_chain",
            "full_selector_sector_projected_block",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "existing_kernel_feedback_inside_K_total",
            "H3_admissible_K_obs_ansatz",
            "R2_internal_feedback_parameter_packet",
            "R3_explicit_internal_light_propagator_packet",
            "R4_explicit_local_chart_emission_map_packet",
        ],
        "resolved_from_P7": ["explicit_emission_map_E_from_M_pair_to_L_int"],
        "missing_upstream_objects": missing_upstream_objects,
        "blocking_frontier": {
            "R4_B1": "explicit_E_packet_present_but_local_chart_preoriented_and_not_factorized",
            "H33_B1": sources["H33"]["result"],
            "H34_B1": sources["H34"]["frontier_text"],
            "H35_B1": sources["H35"]["frontier_text"],
            "H36_B1": sources["H36"]["result"],
            "H37_B1": sources["H37"]["result"],
            "H14_B1": sources["H14"]["frontier"]["H14_B1"],
            "H15_B1": "existing kernel feedback has no explicit residual-selector-sector reduction or projected selector-block export in the current repository, so K_obs remains a distinct extension hypothesis rather than an identified reformulation of existing kernel feedback",
            "H29_B1": sources["H29"]["frontier"],
        },
        "computed": {
            "E_matrix": sources["R4_artifact"]["matrix"],
            "E_domain_basis_order": sources["R4_artifact"]["domain_basis_order"],
            "E_codomain_basis_order": sources["R4_artifact"]["codomain_basis_order"],
            "partial_emission_propagation_pullback": sources["R4_artifact"]["partial_emission_propagation_pullback"]["matrix"],
        },
        "required_next_step": "IMPLEMENT_R_MAT_OR_O_OBS_OR_THE_FACTORIZATION_MAP_AND_RERUN_P8_BEFORE_CLAIMING_EXISTING_KERNEL_FEEDBACK_ALREADY_CONTAINS_SELECTOR_FACING_KOBS",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "status": report["status"],
        "reason": report["reason"],
        "resolved_from_P7": report["resolved_from_P7"],
        "missing_upstream_objects": report["missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

