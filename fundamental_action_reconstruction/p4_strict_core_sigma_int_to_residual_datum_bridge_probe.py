#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p4_strict_core_sigma_int_to_residual_datum_bridge_probe.json"
OUT_SUMMARY = GENERATED / "p4_strict_core_sigma_int_to_residual_datum_bridge_probe_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def blocker_present(items: list[str], needle: str) -> bool:
    return needle in items


def main() -> None:
    sources = {
        "B4": load_json("fundamental_action_reconstruction/generated/b4_minimal_sigma_int_candidate_summary.json"),
        "B5": load_json("fundamental_action_reconstruction/generated/b5_sigma_int_local_stability_audit_summary.json"),
        "B6": load_json("fundamental_action_reconstruction/generated/b6_sigma_to_selector_factorized_bridge_summary.json"),
        "B7": load_json("fundamental_action_reconstruction/generated/b7_factorized_selector_mode_scaffold_compatibility_audit_summary.json"),
        "B8": load_json("fundamental_action_reconstruction/generated/b8_selector_track_anti_overclaim_audit_summary.json"),
        "C37": load_json("fundamental_action_reconstruction/generated/c37_residual_orientation_datum_internalization_candidate_audit_summary.json"),
        "C40": load_json("fundamental_action_reconstruction/generated/c40_minimal_field_list_audit_summary.json"),
        "C46": load_json("fundamental_action_reconstruction/generated/c46_minimal_template_file_creation_audit_summary.json"),
        "T2": load_json("fundamental_action_reconstruction/generated/t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"),
        "AX3": load_json("fundamental_action_reconstruction/generated/ax3_axiom_lane_sigma_int_residual_datum_bridge_instance_summary.json"),
    }

    b8_blockers = sources["B8"]["residual_blockers"]
    c37_result = sources["C37"]["result"]
    ax3_result = sources["AX3"]["result"]

    route_checks = [
        {
            "id": "b4_sigma_candidate_exists",
            "pass": sources["B4"]["b4"]["candidate_name"] == "sigma_int_candidate",
            "expected": "sigma_int_candidate",
            "actual": sources["B4"]["b4"]["candidate_name"],
            "meaning": "the route starts from a named internal candidate object",
        },
        {
            "id": "b8_strict_derivation_blocker_active",
            "pass": blocker_present(b8_blockers, "no_strict_derivation_of_sigma_int_candidate"),
            "expected": True,
            "actual": blocker_present(b8_blockers, "no_strict_derivation_of_sigma_int_candidate"),
            "meaning": "the current route still lacks strict derivation of sigma_int_candidate",
        },
        {
            "id": "b5_gauge_quotient_safety_open",
            "pass": sources["B5"]["b5"]["findings"][2]["status"] == "open",
            "expected": "open",
            "actual": sources["B5"]["b5"]["findings"][2]["status"],
            "meaning": "full gauge-quotient safety is still open",
        },
        {
            "id": "b6_candidate_fit_present",
            "pass": sources["B6"]["findings"]["sigma_fits_residual_z2_orientation_slot"]["status"] == "supported_candidate_fit",
            "expected": "supported_candidate_fit",
            "actual": sources["B6"]["findings"]["sigma_fits_residual_z2_orientation_slot"]["status"],
            "meaning": "sigma_int_candidate fits the residual Z2 slot only as candidate-fit",
        },
        {
            "id": "b7_overlay_only",
            "pass": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"] == "partial_control_route_only",
            "expected": "partial_control_route_only",
            "actual": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"],
            "meaning": "compatibility remains overlay-only rather than strict-core bridge discharge",
        },
        {
            "id": "c37_candidate_internalization_only",
            "pass": c37_result["candidate_internalization_present"] == "yes_candidate_fit",
            "expected": "yes_candidate_fit",
            "actual": c37_result["candidate_internalization_present"],
            "meaning": "internalization is present only at candidate-fit level",
        },
        {
            "id": "c37_strict_core_export_absent",
            "pass": c37_result["strict_core_residual_orientation_datum_export"] == "not_shown",
            "expected": "not_shown",
            "actual": c37_result["strict_core_residual_orientation_datum_export"],
            "meaning": "strict core still does not export the residual orientation datum",
        },
        {
            "id": "t2_map_absent",
            "pass": sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"] is False,
            "expected": False,
            "actual": sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"],
            "meaning": "the strict-core equivalence/export map is absent",
        },
        {
            "id": "c40_semantic_target_field_present",
            "pass": bool(sources["C40"]["field_list"]["target_slot_or_target_datum"]),
            "expected": True,
            "actual": bool(sources["C40"]["field_list"]["target_slot_or_target_datum"]),
            "meaning": "a semantic target field exists in the carrier grammar, but that is not yet a strict-core export",
        },
        {
            "id": "c46_acceptance_carrier_present",
            "pass": bool(sources["C46"]["created_file"]["exists_after_step"]),
            "expected": True,
            "actual": bool(sources["C46"]["created_file"]["exists_after_step"]),
            "meaning": "a persisted acceptance carrier exists, but it is not itself a bridge theorem or export",
        },
        {
            "id": "ax3_axiom_lane_bridge_only",
            "pass": ax3_result["sigma_int_bridge_instance_available"] == "yes_axiom_lane_only",
            "expected": "yes_axiom_lane_only",
            "actual": ax3_result["sigma_int_bridge_instance_available"],
            "meaning": "an explicit bridge instance exists only on the axiom-augmented lane",
        },
        {
            "id": "ax3_does_not_change_strict_core",
            "pass": ax3_result["strict_core_changed"] is False,
            "expected": False,
            "actual": ax3_result["strict_core_changed"],
            "meaning": "the axiom-lane bridge instance does not promote the bridge into strict core",
        },
    ]

    route_state = {
        "sigma_int_candidate_exists": sources["B4"]["b4"]["candidate_name"] == "sigma_int_candidate",
        "strict_core_source_object_present": not blocker_present(b8_blockers, "no_strict_derivation_of_sigma_int_candidate"),
        "theorem_level_gauge_quotient_safety_present": sources["B5"]["b5"]["findings"][2]["status"] != "open",
        "candidate_fit_present": sources["B6"]["findings"]["sigma_fits_residual_z2_orientation_slot"]["status"] == "supported_candidate_fit",
        "semantic_target_field_present": bool(sources["C40"]["field_list"]["target_slot_or_target_datum"]),
        "acceptance_artifact_carrier_present": bool(sources["C46"]["created_file"]["exists_after_step"]),
        "axiom_lane_bridge_instance_present": ax3_result["sigma_int_bridge_instance_available"] == "yes_axiom_lane_only",
        "strict_core_target_slot_export_present": c37_result["strict_core_residual_orientation_datum_export"] != "not_shown",
        "strict_core_equivalence_or_export_map_present": bool(sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"]),
        "selector_track_identification_beyond_overlay_present": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"] != "partial_control_route_only",
        "strict_core_bridge_discharged": False,
    }

    missing_upstream_objects: list[str] = []
    if not route_state["strict_core_source_object_present"]:
        missing_upstream_objects.append("strict_derivation_or_source_object_upgrade_for_sigma_int_candidate")
    if not route_state["theorem_level_gauge_quotient_safety_present"]:
        missing_upstream_objects.append("theorem_level_gauge_quotient_safety_for_sigma_int_candidate")
    if not route_state["strict_core_target_slot_export_present"]:
        missing_upstream_objects.append("strict_core_target_slot_export_for_residual_orientation_datum")
    if not route_state["strict_core_equivalence_or_export_map_present"]:
        missing_upstream_objects.append(
            "strict_core_equivalence_or_export_map_sigma_int_candidate_to_residual_orientation_datum"
        )
    if not route_state["selector_track_identification_beyond_overlay_present"]:
        missing_upstream_objects.append("selector_track_identification_beyond_overlay_only_for_sigma_int_bridge")

    report = {
        "stage": "P4",
        "goal": "compute_or_fail_strict_core_sigma_int_to_residual_orientation_datum",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE",
        "reason": "the current route stops at candidate-fit, acceptance carrier, and axiom-lane-only bridge witness before a strict-core residual-datum export or bridge map exists",
        "lane": "strict_core_sigma_int_residual_datum_route",
        "route_under_test": [
            "sigma_int_candidate",
            "residual_orientation_datum",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "sigma_int_candidate",
            "residual_Z2_candidate_fit",
            "semantic_target_field_in_acceptance_carrier",
            "persisted_acceptance_artifact_carrier",
            "axiom_lane_sigma_int_bridge_instance",
        ],
        "missing_upstream_objects": missing_upstream_objects,
        "blocking_frontier": {
            "B8": "no_strict_derivation_of_sigma_int_candidate / no_theorem_level_gauge_quotient_safety",
            "C37_B1": sources["C37"]["residual_blockers"]["C37_B1"],
            "T2_B1": sources["T2"]["frontier_after_T2"]["T2_B1"],
            "AX3_result": "sigma_int_bridge_instance_is_materialized_on_axiom_lane_only_and_does_not_change_strict_core",
        },
        "computed": {},
        "required_next_step": "IMPLEMENT_ONE_MISSING_STRICT_CORE_RESIDUAL_DATUM_BRIDGE_OBJECT_AND_RERUN_P4_BEFORE_CLAIMING_STRICT_CORE_INTERNALIZATION",
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
