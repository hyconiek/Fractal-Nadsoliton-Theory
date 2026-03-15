#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p2_strict_core_sigma_int_to_a1_pair1_probe.json"
OUT_SUMMARY = GENERATED / "p2_strict_core_sigma_int_to_a1_pair1_probe_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "QW2191": load_json(
            "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
        ),
        "B2": load_json("fundamental_action_reconstruction/generated/b2_internal_orientation_datum_source_audit_summary.json"),
        "B4": load_json("fundamental_action_reconstruction/generated/b4_minimal_sigma_int_candidate_summary.json"),
        "B6": load_json("fundamental_action_reconstruction/generated/b6_sigma_to_selector_factorized_bridge_summary.json"),
        "B7": load_json("fundamental_action_reconstruction/generated/b7_factorized_selector_mode_scaffold_compatibility_audit_summary.json"),
        "T2": load_json("fundamental_action_reconstruction/generated/t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"),
        "C35": load_json("fundamental_action_reconstruction/generated/c35_actual_phase_source_branch_audit_summary.json"),
        "C37": load_json("fundamental_action_reconstruction/generated/c37_residual_orientation_datum_internalization_candidate_audit_summary.json"),
        "C46": load_json("fundamental_action_reconstruction/generated/c46_minimal_template_file_creation_audit_summary.json"),
        "C47": load_json("fundamental_action_reconstruction/generated/c47_basis_level_orientation_slice_candidate_audit_summary.json"),
        "C48": load_json("fundamental_action_reconstruction/generated/c48_minimal_actual_basis_pair_export_skeleton_audit_summary.json"),
        "C49": load_json("fundamental_action_reconstruction/generated/c49_conditional_populated_instance_schema_audit_summary.json"),
    }

    route_checks = [
        {
            "id": "qw2191_extra_symmetry_breaking_required",
            "pass": bool(sources["QW2191"]["flags"]["obstruction_requires_explicit_symmetry_breaking_postulate"]),
            "expected": True,
            "actual": bool(sources["QW2191"]["flags"]["obstruction_requires_explicit_symmetry_breaking_postulate"]),
            "meaning": "kernel alone requires extra symmetry breaking before uniqueness can close",
        },
        {
            "id": "b4_sigma_candidate_exists",
            "pass": sources["B4"]["b4"]["candidate_name"] == "sigma_int_candidate",
            "expected": "sigma_int_candidate",
            "actual": sources["B4"]["b4"]["candidate_name"],
            "meaning": "a minimal internal datum candidate exists",
        },
        {
            "id": "b2_internal_orientation_datum_axis_only_present",
            "pass": sources["B2"]["b2"]["strict_internal_selector_derivations_found"] > 0,
            "expected": ">0",
            "actual": sources["B2"]["b2"]["strict_internal_selector_derivations_found"],
            "meaning": "an axis-only internal orientation datum is exported on at least one declared lane (residual Z2 remains explicit)",
        },
        {
            "id": "b6_residual_z2_fit_present",
            "pass": sources["B6"]["findings"]["sigma_fits_residual_z2_orientation_slot"]["status"] == "supported_candidate_fit",
            "expected": "supported_candidate_fit",
            "actual": sources["B6"]["findings"]["sigma_fits_residual_z2_orientation_slot"]["status"],
            "meaning": "sigma_int_candidate fits the residual Z2 slot only at candidate-fit level",
        },
        {
            "id": "b7_overlay_only",
            "pass": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"] == "partial_control_route_only",
            "expected": "partial_control_route_only",
            "actual": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"],
            "meaning": "the bridge remains overlay/control-route only and not strict-core discharge",
        },
        {
            "id": "t2_strict_core_export_map_present",
            "pass": bool(sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"]),
            "expected": True,
            "actual": bool(sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"]),
            "meaning": "a strict-core sigma-int -> residual export-map object exists in declared R1 scope (no implied selector closure)",
        },
        {
            "id": "c35_strict_core_actual_phase_source_present",
            "pass": sources["C35"]["result"]["strict_core_actual_phase_source"] != "not_shown",
            "expected": True,
            "actual": sources["C35"]["result"]["strict_core_actual_phase_source"] != "not_shown",
            "meaning": "strict core exports an actual theta-source on declared pair1/pair2 lanes (residual Z2 remains explicit)",
        },
        {
            "id": "c37_candidate_internalization_only",
            "pass": sources["C37"]["result"]["candidate_internalization_present"] == "yes_candidate_fit",
            "expected": "yes_candidate_fit",
            "actual": sources["C37"]["result"]["candidate_internalization_present"],
            "meaning": "internalization is still only candidate-fit, not strict equivalence",
        },
        {
            "id": "c46_acceptance_carrier_exists",
            "pass": bool(sources["C46"]["created_file"]["exists_after_step"]),
            "expected": True,
            "actual": bool(sources["C46"]["created_file"]["exists_after_step"]),
            "meaning": "an acceptance-artifact carrier exists, but that is not yet a bridge theorem/export",
        },
        {
            "id": "c47_orientation_slice_class_present",
            "pass": sources["C47"]["findings"]["class_level_basis_candidate"] == "present_partial",
            "expected": "present_partial",
            "actual": sources["C47"]["findings"]["class_level_basis_candidate"],
            "meaning": "a class-level orientation-slice candidate exists",
        },
        {
            "id": "c48_basis_pair_export_skeleton_present",
            "pass": sources["C48"]["findings"]["minimal_export_skeleton"] == "present_partial",
            "expected": "present_partial",
            "actual": sources["C48"]["findings"]["minimal_export_skeleton"],
            "meaning": "a minimal basis-pair export skeleton exists",
        },
        {
            "id": "c49_actual_theta_supply_present",
            "pass": sources["C49"]["findings"]["actual_theta_supply"] != "not_shown",
            "expected": True,
            "actual": sources["C49"]["findings"]["actual_theta_supply"] != "not_shown",
            "meaning": "the basis-pair conditional schema has strict-core theta supply in declared scope",
        },
    ]

    route_state = {
        "sigma_int_candidate_exists": sources["B4"]["b4"]["candidate_name"] == "sigma_int_candidate",
        "strict_core_source_object_present": sources["B2"]["b2"]["strict_internal_selector_derivations_found"] > 0,
        "residual_z2_candidate_fit_present": sources["B6"]["findings"]["sigma_fits_residual_z2_orientation_slot"]["status"] == "supported_candidate_fit",
        "strict_core_equivalence_or_export_map_present": bool(sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"]),
        "strict_core_actual_theta_supply_present": sources["C35"]["result"]["strict_core_actual_phase_source"] != "not_shown",
        "acceptance_artifact_carrier_exists": bool(sources["C46"]["created_file"]["exists_after_step"]),
        "orientation_slice_candidate_class_present": sources["C47"]["findings"]["class_level_basis_candidate"] == "present_partial",
        "basis_pair_export_skeleton_present": sources["C48"]["findings"]["minimal_export_skeleton"] == "present_partial",
        "actual_basis_pair_export_present": sources["C48"]["findings"]["actual_basis_pair_export"] != "not_shown",
        "conditional_basis_pair_population_schema_present": sources["C49"]["findings"]["conditional_populated_instance_schema"] == "present_partial",
        "strict_core_route_reaches_operator_stage": False,
    }

    missing_upstream_objects: list[str] = []
    if not route_state["strict_core_source_object_present"]:
        missing_upstream_objects.append(
            "strict_core_source_object_upgrading_sigma_int_candidate_beyond_candidate_only_status"
        )
    if not route_state["strict_core_equivalence_or_export_map_present"]:
        missing_upstream_objects.append(
            "strict_core_equivalence_or_export_map_sigma_int_candidate_to_residual_orientation_datum"
        )
    if not route_state["strict_core_actual_theta_supply_present"]:
        missing_upstream_objects.append(
            "strict_core_actual_theta_1_theta_2_supply_for_current_pair_frames"
        )
    if not route_state["actual_basis_pair_export_present"]:
        missing_upstream_objects.append(
            "populated_actual_basis_pair_export_u_1_u_2_for_current_pair_frames"
        )
    missing_upstream_objects.append(
        "strict_core_operator_level_bridge_from_materialized_orientation_slice_to_A_1_pair1"
    )

    report = {
        "stage": "P2",
        "goal": "compute_or_fail_strict_core_sigma_int_to_A1_pair1",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_ROUTE",
        "reason": "the route reaches strict-core theta supply and basis-pair materialization in declared scope, but no strict-core operator-level export to A_1(pair1) is yet exported",
        "lane": "strict_core_candidate_route",
        "route_under_test": [
            "sigma_int_candidate",
            "residual_orientation_datum",
            "theta_1_theta_2",
            "u_1_u_2_and_orientation_slice",
            "A_1(pair1)",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "sigma_int_candidate",
            "residual_Z2_candidate_fit",
            "acceptance_artifact_carrier",
            "orientation_slice_candidate_class",
            "basis_pair_export_skeleton",
            "conditional_basis_pair_population_schema",
            "strict_core_theta_supply_in_declared_scope",
            "strict_core_target_slot_population_in_declared_scope",
        ],
        "missing_upstream_objects": missing_upstream_objects,
        "blocking_frontier": {
            "B2_reduced_blocker": sources["B2"]["b2"]["reduced_blocker"],
            "T2_B1": sources["T2"]["frontier_after_T2"]["T2_B1"],
            "C35_B1": sources["C35"]["residual_blockers"]["C35_B1"],
            "C49_B1": sources["C49"]["frontier_after_C49"]["C49_B1"],
            "C32_B2": sources["C49"]["frontier_after_C49"]["C32_B2"],
        },
        "computed": {},
        "required_next_step": "IMPLEMENT_ONE_OF_THE_MISSING_STRICT_CORE_OBJECTS_AND_RERUN_P2_BEFORE_CLAIMING_A1_PAIR1_STRICT_CORE_REACHABILITY",
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
