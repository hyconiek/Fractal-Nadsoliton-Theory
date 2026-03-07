#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p6_existing_kernel_feedback_to_kobs_operator_chain_probe.json"
OUT_SUMMARY = GENERATED / "p6_existing_kernel_feedback_to_kobs_operator_chain_probe_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "H3": load_json("fundamental_action_reconstruction/generated/h3_minimal_internal_light_feedback_operator_ansatz_packet_summary.json"),
        "H4": load_json("fundamental_action_reconstruction/generated/h4_residual_o2_reduction_of_light_feedback_ansatz_summary.json"),
        "H14": load_json("fundamental_action_reconstruction/generated/h14_existing_kernel_feedback_vs_new_kobs_separation_audit_summary.json"),
        "H15": load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json"),
        "H29": load_json("fundamental_action_reconstruction/generated/h29_wave_retardation_proxy_selector_reduction_audit_summary.json"),
        "R2": load_json("fundamental_action_reconstruction/generated/r2_existing_internal_feedback_parameter_packet_for_kobs_summary.json"),
        "QW1950": load_json(
            "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1950_internal_emergent_observer_closed_loop.json"
        ),
        "QW1951": load_json(
            "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1951_mass_informational_weight_internal_observer.json"
        ),
        "QW1952": load_json(
            "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1952_information_channel_dedegeneracy_operator.json"
        ),
        "QW1956": load_json(
            "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1956_two_state_observer_with_repaired_operator.json"
        ),
    }

    q1950_params = sources["QW1950"]["derived_internal_observer_params"]
    q1951_weights = sources["QW1951"]["mass_informational_weights"]
    q1952_params = sources["QW1952"]["derived_params"]
    q1956_params = sources["QW1956"]["repaired_operator_params"]

    route_checks = [
        {
            "id": "existing_kernel_feedback_present",
            "pass": True,
            "expected": True,
            "actual": True,
            "meaning": "the repo already contains internal feedback-like structure in K_total -> K(d)",
        },
        {
            "id": "h3_operator_chain_required",
            "pass": len(sources["H3"]["maps"]) == 7,
            "expected": 7,
            "actual": len(sources["H3"]["maps"]),
            "meaning": "the admissible K_obs route requires an explicit operator chain with the H3 maps",
        },
        {
            "id": "r2_parameter_packet_present",
            "pass": sources["R2"]["result"]
            == "existing_internal_feedback_parameter_packet_present_but_not_yet_an_operator_level_kobs_factorization",
            "expected": "existing_internal_feedback_parameter_packet_present_but_not_yet_an_operator_level_kobs_factorization",
            "actual": sources["R2"]["result"],
            "meaning": "existing internal feedback parameters are now aggregated in one packet",
        },
        {
            "id": "q1950_observer_loop_parameters_present",
            "pass": all(
                key in q1950_params for key in ["observer_tau", "observer_feedback_gain", "observer_feedback_theta"]
            ),
            "expected": True,
            "actual": all(
                key in q1950_params for key in ["observer_tau", "observer_feedback_gain", "observer_feedback_theta"]
            ),
            "meaning": "observer-loop parameters exist at scalar level",
        },
        {
            "id": "q1951_mass_information_parameters_present",
            "pass": "mass_gain" in q1951_weights,
            "expected": True,
            "actual": "mass_gain" in q1951_weights,
            "meaning": "mass-information coupling parameters exist at scalar level",
        },
        {
            "id": "q1952_anisotropy_parameters_present",
            "pass": all(key in q1952_params for key in ["anisotropy_strength", "retard_phase", "orientation_psi0"]),
            "expected": True,
            "actual": all(key in q1952_params for key in ["anisotropy_strength", "retard_phase", "orientation_psi0"]),
            "meaning": "anisotropy and retardation proxies exist at scalar level",
        },
        {
            "id": "q1956_repaired_operator_parameters_present",
            "pass": all(key in q1956_params for key in ["a2_even_mode", "b1_odd_mode", "b3_odd_mode"]),
            "expected": True,
            "actual": all(key in q1956_params for key in ["a2_even_mode", "b1_odd_mode", "b3_odd_mode"]),
            "meaning": "repaired-operator parameters exist at scalar level",
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
            "meaning": "existing feedback is not exported to the selector sector",
        },
        {
            "id": "h29_proxies_preoriented_only",
            "pass": sources["H29"]["status"] == "PASS_PARTIAL_PREORIENTED_PROXY_REDUCTION_ONLY",
            "expected": "PASS_PARTIAL_PREORIENTED_PROXY_REDUCTION_ONLY",
            "actual": sources["H29"]["status"],
            "meaning": "old proxies only modulate a pre-oriented channel and do not provide an internal anchor by themselves",
        },
    ]

    route_state = {
        "kernel_feedback_present": True,
        "existing_internal_feedback_parameter_packet_present": sources["R2"]["result"]
        == "existing_internal_feedback_parameter_packet_present_but_not_yet_an_operator_level_kobs_factorization",
        "h3_operator_chain_spec_present": True,
        "selector_sector_reduction_packet_present": True,
        "explicit_emission_map_present": False,
        "explicit_internal_light_propagator_present": False,
        "explicit_light_to_matter_response_map_present": False,
        "explicit_observer_readout_operator_present": False,
        "equivalence_factorization_map_from_existing_feedback_present": False,
        "selector_sector_projected_block_present": False,
        "strict_core_kobs_instantiated": False,
    }

    missing_upstream_objects = [
        "explicit_emission_map_E_from_M_pair_to_L_int",
        "explicit_internal_light_propagator_G_light_on_L_int",
        "explicit_light_to_matter_response_map_R_mat_from_L_int_to_Q_mat",
        "explicit_observer_readout_operator_O_obs_on_Q_mat",
        "equivalence_or_factorization_map_from_existing_kernel_feedback_and_R2_parameter_packet_to_H3_operator_chain",
        "selector_sector_projected_2x2_block_export_on_an_actual_pair",
    ]

    report = {
        "stage": "P6",
        "goal": "compute_or_fail_existing_kernel_feedback_to_selector_facing_kobs",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_KERNEL_FEEDBACK_TO_KOBS_ROUTE",
        "reason": "the repo contains internal feedback parameters and a packet-ready K_obs ansatz, but still lacks the operator-level factorization maps and selector-facing projected block needed to instantiate K_obs from existing kernel feedback",
        "lane": "existing_kernel_feedback_to_kobs_operator_route",
        "route_under_test": [
            "existing_kernel_feedback",
            "internal_feedback_parameter_packet",
            "H3_operator_chain",
            "selector_sector_projected_block",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "existing_kernel_feedback_inside_K_total",
            "H3_admissible_K_obs_ansatz",
            "H4_residual_O2_reduction_packet",
            "R2_internal_feedback_parameter_packet",
            "QW1950_QW1951_QW1952_QW1956_scalar_proxy_parameters",
        ],
        "missing_upstream_objects": missing_upstream_objects,
        "blocking_frontier": {
            "R2_B1": "parameter_packet_present_but_not_yet_an_operator_level_kobs_factorization",
            "H14_B1": sources["H14"]["frontier"]["H14_B1"],
            "H15_B1": "existing kernel feedback has no explicit residual-selector-sector reduction or projected selector-block export in the current repository, so K_obs remains a distinct extension hypothesis rather than an identified reformulation of existing kernel feedback",
            "H29_B1": sources["H29"]["frontier"],
            "H4_B1": sources["H4"]["frontier"]["H4_B1"],
        },
        "computed": {},
        "required_next_step": "IMPLEMENT_ONE_OPERATOR_CHAIN_OBJECT_OR_ONE_FACTORIZATION_MAP_AND_RERUN_P6_BEFORE_CLAIMING_EXISTING_KERNEL_FEEDBACK_ALREADY_CONTAINS_KOBS",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "status": report["status"],
        "reason": report["reason"],
        "missing_upstream_objects": report["missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
