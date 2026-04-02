#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1117 = GENERATED / "p1117_current_strict_t173_t176_f973_f960_no_exported_live_nonsamelane_provider_audit_probe_summary.json"
IN_F974 = GENERATED / "f974_current_strict_t173_t176_f973_f960_no_exported_live_nonsamelane_provider_packet_summary.json"
IN_P721 = GENERATED / "p721_current_strict_t176_source_topology_basis_free_qw2191_safe_provider_nonupgrade_audit_probe_summary.json"
IN_P722 = GENERATED / "p722_current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_audit_probe_summary.json"
IN_P723 = GENERATED / "p723_current_strict_t178_source_topology_to_atlas_chart_seed_selection_bridge_nonexport_audit_probe_summary.json"
IN_P724 = GENERATED / "p724_current_strict_t178_positive_source_polarity_atlas_entry_corridor_reduction_audit_probe_summary.json"
IN_P725 = GENERATED / "p725_current_strict_t179_positive_corridor_odd_even_lane_selection_bridge_nonexport_audit_probe_summary.json"
IN_P726 = GENERATED / "p726_current_strict_t180_positive_corridor_outer_interior_chart_selection_bridge_nonexport_audit_probe_summary.json"
IN_P1059 = GENERATED / "p1059_current_strict_t173_t176_source_side_input_leg_same_lane_stagnation_and_stop_audit_probe_summary.json"

OUT_JSON = GENERATED / "p1118_current_strict_t173_t176_f974_f960_t178_chart_seed_candidate_admission_probe.json"
OUT_SUMMARY = GENERATED / "p1118_current_strict_t173_t176_f974_f960_t178_chart_seed_candidate_admission_probe_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
T178_TARGET = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
ACTIVE_CONTRACT = "attack_existing_t178_chart_seed_candidate_beneath_existing_f960_bridge_target"
CANDIDATE_GRADE = "current_best_new_narrow_provider_class_seed_only"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1117, IN_F974, IN_P721, IN_P722, IN_P723, IN_P724, IN_P725, IN_P726, IN_P1059]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1118",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1117 = load_json(IN_P1117)
    f974 = load_json(IN_F974)
    p721 = load_json(IN_P721)
    p722 = load_json(IN_P722)
    p723 = load_json(IN_P723)
    p724 = load_json(IN_P724)
    p725 = load_json(IN_P725)
    p726 = load_json(IN_P726)
    p1059 = load_json(IN_P1059)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append({
            "id": check_id,
            "actual": actual,
            "expected": expected,
            "pass": passed,
            "meaning": meaning,
        })
        if not passed:
            blocking.append(check_id)

    no_live_candidate_frozen = (
        p1117.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F973_F960_NO_ALREADY_EXPORTED_LIVE_NON_SAME_LANE_INVERSION_SENSITIVE_SOURCE_SIDE_PROVIDER_CANDIDATE_AUDITED"
        and p1117.get("already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present") is False
        and p1117.get("active_missing_bridge") == ACTIVE_BRIDGE
        and f974.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F973_F960_NO_ALREADY_EXPORTED_LIVE_NON_SAME_LANE_INVERSION_SENSITIVE_SOURCE_SIDE_PROVIDER_CANDIDATE_PACKET_EXPORTED"
        and f974.get("packet_exported_on_current_repo_state") is True
    )

    source_topology_lane_chart_blind_gap_frozen = (
        p721.get("status") == "PASS_SOURCE_TOPOLOGY_BASIS_FREE_QW2191_SAFE_PROVIDER_NONUPGRADE_AUDITED"
        and p721.get("source_topology_physically_interpretable_strict_ingredients_present") is True
        and p721.get("source_topology_lane_upgrades_to_exact_t176_provider") is False
        and p721.get("preferred_next_direction") == "chart_sensitive_transported_flux_or_current_like_provider_class_over_further_static_basis_free_or_output_axis_classes"
        and p722.get("status") == "PASS_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_NONEXPORT_AUDITED"
        and p722.get("current_source_topology_lane_is_physics_facing_but_chart_blind") is True
        and p722.get("t177_target_exported_on_current_repo_state") is False
        and p722.get("next_honest_move") == "chart_sensitive_transported_flux_current_like_section_bridge"
    )

    existing_t178_nonexport_target_is_exact_next_seed_gap = (
        p723.get("status") == "PASS_SOURCE_TO_ATLAS_CHART_SEED_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p723.get("t178_target_name") == T178_TARGET
        and p723.get("t178_target_exported_on_current_repo_state") is False
        and p723.get("current_source_topology_lane_supplies_sign_flow_and_selector_polarity_but_not_chart_seed_selection") is True
        and p723.get("next_honest_move") == "source_topology_to_atlas_chart_seed_selection_bridge"
    )

    chart_seed_gap_still_unresolved = (
        p724.get("status") == "PARTIAL_POSITIVE_SOURCE_POLARITY_REDUCES_ATLAS_ENTRY_CORRIDOR_ONLY"
        and p724.get("source_positive_polarity_available") is True
        and p724.get("unique_chart_seed_selected") is False
        and p725.get("status") == "PASS_POSITIVE_CORRIDOR_ODD_EVEN_LANE_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p725.get("t179_target_exported_on_current_repo_state") is False
        and p726.get("status") == "PASS_POSITIVE_CORRIDOR_OUTER_INTERIOR_CHART_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p726.get("t180_target_exported_on_current_repo_state") is False
    )

    stopped_source_side_input_leg_same_lane_not_reactivated = (
        p1059.get("status") == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_SAME_LANE_STAGNATION_AND_STOP_AUDITED"
        and p1059.get("same_lane_stagnation_boundary_reached") is True
        and p1059.get("further_same_lane_descent_is_not_honest_primary_move") is True
        and p1059.get("restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route") is True
    )

    candidate_admissible_against_existing_f960_bridge_target = (
        no_live_candidate_frozen
        and source_topology_lane_chart_blind_gap_frozen
        and existing_t178_nonexport_target_is_exact_next_seed_gap
        and chart_seed_gap_still_unresolved
        and stopped_source_side_input_leg_same_lane_not_reactivated
    )

    add_check("no_live_candidate_frozen", no_live_candidate_frozen, True, "P1117/F974 already freeze the absence of any already-live non-same-lane inversion-sensitive source-side provider candidate beneath F960.")
    add_check("source_topology_lane_chart_blind_gap_frozen", source_topology_lane_chart_blind_gap_frozen, True, "P721/N717/P722 already freeze that the source-topology lane is physics-facing but chart-blind and points toward a chart-sensitive transported-section family.")
    add_check("existing_t178_nonexport_target_is_exact_next_seed_gap", existing_t178_nonexport_target_is_exact_next_seed_gap, True, "P723/N719 already freeze SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1 as the sharp exact nonexport chart-seed gap.")
    add_check("chart_seed_gap_still_unresolved", chart_seed_gap_still_unresolved, True, "P724/P725/P726 already show that the chart-seed gap remains unresolved and still lacks a unique atlas seed selection.")
    add_check("stopped_source_side_input_leg_same_lane_not_reactivated", stopped_source_side_input_leg_same_lane_not_reactivated, True, "P1059 already freezes that this admission cannot be read as restarting the stopped source-side input-leg same-lane family.")
    add_check("candidate_admissible_against_existing_f960_bridge_target", candidate_admissible_against_existing_f960_bridge_target, True, "Therefore the current best new narrow provider-class seed candidate beneath F960 may now be admitted as the existing T178 chart-seed bridge target only.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F974_F960_EXISTING_T178_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_ADMITTED"
        if not blocking else
        "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F974_F960_EXISTING_T178_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_ADMISSION"
    )

    artifact = {
        "stage": "P1118",
        "status": status,
        "as_of": AS_OF,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "current_best_new_narrow_provider_class_seed_candidate": T178_TARGET,
        "candidate_admissible_against_existing_f960_bridge_target": candidate_admissible_against_existing_f960_bridge_target,
        "candidate_grade": CANDIDATE_GRADE,
        "candidate_is_already_exported_live_provider": False,
        "candidate_counts_as_actual_t178_export": False,
        "counts_as_lawful_supplier": False,
        "counts_as_solution": False,
        "counts_as_strict_physical_orientation_datum": False,
        "active_contract": ACTIVE_CONTRACT,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "current_best_new_narrow_provider_class_seed_candidate": artifact["current_best_new_narrow_provider_class_seed_candidate"],
        "candidate_admissible_against_existing_f960_bridge_target": artifact["candidate_admissible_against_existing_f960_bridge_target"],
        "candidate_grade": artifact["candidate_grade"],
        "active_contract": artifact["active_contract"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
