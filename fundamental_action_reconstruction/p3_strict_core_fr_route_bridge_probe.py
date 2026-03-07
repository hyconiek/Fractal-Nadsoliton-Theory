#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p3_strict_core_fr_route_bridge_probe.json"
OUT_SUMMARY = GENERATED / "p3_strict_core_fr_route_bridge_probe_summary.json"


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
        "T2": load_json("fundamental_action_reconstruction/generated/t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"),
        "C35": load_json("fundamental_action_reconstruction/generated/c35_actual_phase_source_branch_audit_summary.json"),
        "C37": load_json("fundamental_action_reconstruction/generated/c37_residual_orientation_datum_internalization_candidate_audit_summary.json"),
        "C38": load_json("fundamental_action_reconstruction/generated/c38_sigma_int_residual_datum_theorem_spec_audit_summary.json"),
        "C46": load_json("fundamental_action_reconstruction/generated/c46_minimal_template_file_creation_audit_summary.json"),
    }

    b8_blockers = sources["B8"]["residual_blockers"]

    route_checks = [
        {
            "id": "b4_sigma_candidate_exists",
            "pass": sources["B4"]["b4"]["candidate_name"] == "sigma_int_candidate",
            "expected": "sigma_int_candidate",
            "actual": sources["B4"]["b4"]["candidate_name"],
            "meaning": "the FR route has a canonical candidate datum",
        },
        {
            "id": "b8_no_strict_derivation_of_sigma",
            "pass": blocker_present(b8_blockers, "no_strict_derivation_of_sigma_int_candidate"),
            "expected": True,
            "actual": blocker_present(b8_blockers, "no_strict_derivation_of_sigma_int_candidate"),
            "meaning": "sigma_int_candidate is not yet strict-derived",
        },
        {
            "id": "b5_full_gauge_quotient_safety_open",
            "pass": sources["B5"]["b5"]["findings"][2]["status"] == "open",
            "expected": "open",
            "actual": sources["B5"]["b5"]["findings"][2]["status"],
            "meaning": "full gauge-quotient safety remains open",
        },
        {
            "id": "b6_sigma_alone_to_theta_not_shown",
            "pass": sources["B6"]["findings"]["sigma_alone_selects_theta"]["status"] == "not_shown",
            "expected": "not_shown",
            "actual": sources["B6"]["findings"]["sigma_alone_selects_theta"]["status"],
            "meaning": "sigma alone does not yet derive theta",
        },
        {
            "id": "b6_factorized_bridge_only",
            "pass": sources["B6"]["findings"]["factorized_bridge"]["status"] == "candidate_control_bridge_identified",
            "expected": "candidate_control_bridge_identified",
            "actual": sources["B6"]["findings"]["factorized_bridge"]["status"],
            "meaning": "the bridge exists only as a factorized control route",
        },
        {
            "id": "b7_overlay_compatibility_only",
            "pass": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"] == "partial_control_route_only",
            "expected": "partial_control_route_only",
            "actual": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"],
            "meaning": "compatibility is overlay-only, not strict-core discharge",
        },
        {
            "id": "t2_map_absent",
            "pass": sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"] is False,
            "expected": False,
            "actual": sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"],
            "meaning": "no strict-core bridge map exists",
        },
        {
            "id": "c37_equivalence_bridge_not_shown",
            "pass": sources["C37"]["result"]["strict_core_equivalence_bridge"] == "not_shown",
            "expected": "not_shown",
            "actual": sources["C37"]["result"]["strict_core_equivalence_bridge"],
            "meaning": "candidate internalization is not strict equivalence",
        },
        {
            "id": "c38_theorem_spec_absent",
            "pass": sources["C38"]["findings"]["strict_core_theorem_spec_present"] is False,
            "expected": False,
            "actual": sources["C38"]["findings"]["strict_core_theorem_spec_present"],
            "meaning": "no strict-core theorem-spec exists for the bridge",
        },
        {
            "id": "c35_actual_theta_source_absent",
            "pass": sources["C35"]["result"]["strict_core_actual_phase_source"] == "not_shown",
            "expected": "not_shown",
            "actual": sources["C35"]["result"]["strict_core_actual_phase_source"],
            "meaning": "the FR route does not yet supply actual theta_1, theta_2",
        },
        {
            "id": "c46_carrier_exists_but_is_insufficient",
            "pass": bool(sources["C46"]["created_file"]["exists_after_step"]),
            "expected": True,
            "actual": bool(sources["C46"]["created_file"]["exists_after_step"]),
            "meaning": "a carrier exists, but carrier presence is not a bridge theorem/export",
        },
    ]

    route_state = {
        "sigma_int_candidate_exists": True,
        "strict_derivation_of_sigma_int_candidate_present": not blocker_present(
            b8_blockers, "no_strict_derivation_of_sigma_int_candidate"
        ),
        "full_gauge_quotient_safety_present": sources["B5"]["b5"]["findings"][2]["status"] != "open",
        "strict_core_equivalence_or_export_map_present": bool(
            sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"]
        ),
        "standalone_sigma_to_theta_derivation_present": sources["B6"]["findings"]["sigma_alone_selects_theta"]["status"]
        != "not_shown",
        "internal_derivation_of_selector_family_present": not blocker_present(
            b8_blockers, "no_internal_derivation_of_Jab_selector_family"
        ),
        "strict_core_actual_theta_source_present": sources["C35"]["result"]["strict_core_actual_phase_source"]
        != "not_shown",
        "carrier_exists": bool(sources["C46"]["created_file"]["exists_after_step"]),
        "route_reaches_residual_datum": False,
        "route_reaches_theta_source": False,
    }

    missing_upstream_objects: list[str] = []
    if not route_state["strict_derivation_of_sigma_int_candidate_present"]:
        missing_upstream_objects.append(
            "strict_derivation_or_strict_core_source_object_upgrade_for_sigma_int_candidate"
        )
    if not route_state["full_gauge_quotient_safety_present"]:
        missing_upstream_objects.append(
            "theorem_level_gauge_quotient_safety_for_sigma_int_candidate"
        )
    if not route_state["strict_core_equivalence_or_export_map_present"]:
        missing_upstream_objects.append(
            "strict_core_equivalence_or_export_map_sigma_int_candidate_to_residual_orientation_datum"
        )
    if not route_state["standalone_sigma_to_theta_derivation_present"] and not route_state[
        "internal_derivation_of_selector_family_present"
    ]:
        missing_upstream_objects.append(
            "strict_core_sigma_int_to_theta_map_or_internal_derivation_of_Jab_selector_family"
        )
    if not route_state["strict_core_actual_theta_source_present"]:
        missing_upstream_objects.append(
            "strict_core_actual_theta_1_theta_2_source_for_current_pair_frames"
        )

    report = {
        "stage": "P3",
        "goal": "compute_or_fail_strict_core_fr_route_bridge",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_FR_ROUTE",
        "reason": "the FR/topological route remains candidate/control only and does not yet export a strict-core residual-datum bridge or theta source",
        "lane": "strict_core_fr_topological_route",
        "route_under_test": [
            "sigma_int_candidate",
            "residual_orientation_datum",
            "theta_1_theta_2",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "sigma_int_candidate",
            "partial_local_stability_support",
            "factorized_control_bridge",
            "overlay_compatibility",
            "acceptance_artifact_carrier",
        ],
        "missing_upstream_objects": missing_upstream_objects,
        "blocking_frontier": {
            "B8_residual_blockers": b8_blockers,
            "T2_B1": sources["T2"]["frontier_after_T2"]["T2_B1"],
            "C35_B1": sources["C35"]["residual_blockers"]["C35_B1"],
            "C37_B1": sources["C37"]["residual_blockers"]["C37_B1"],
            "C38_B1": sources["C38"]["frontier_after_c38"]["C38_B1"],
        },
        "computed": {},
        "required_next_step": "IMPLEMENT_ONE_MISSING_STRICT_CORE_FR_BRIDGE_OBJECT_AND_RERUN_P3_BEFORE_CLAIMING_INTERNAL_SOURCE_REACHABILITY",
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
