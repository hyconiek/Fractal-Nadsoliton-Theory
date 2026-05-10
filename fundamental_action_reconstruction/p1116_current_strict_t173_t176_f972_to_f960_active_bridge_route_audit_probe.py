#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1115 = GENERATED / "p1115_current_strict_t173_t176_f966_selector_source_subroute_audit_probe_summary.json"
IN_F972 = GENERATED / "f972_current_strict_t173_t176_f966_selector_source_subroute_packet_summary.json"
IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"
IN_P1060 = GENERATED / "p1060_current_strict_t173_t176_post_stop_source_side_input_leg_lane_to_existing_t183_residual_datum_pair12_orbit_direction_frontier_route_decision_audit_probe_summary.json"
IN_P1061 = GENERATED / "p1061_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_existing_export_or_near_miss_candidate_audit_probe_summary.json"
IN_F960 = GENERATED / "f960_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_target_packet_summary.json"

OUT_JSON = GENERATED / "p1116_current_strict_t173_t176_f972_to_f960_active_bridge_route_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1116_current_strict_t173_t176_f972_to_f960_active_bridge_route_audit_probe_summary.json"

PRIMARY_SUBROUTE = "explicit_strict_internal_selector_source_derivation_frontier"
PRIMARY_TARGET = "T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC"
ACTIVE_FRONTIER = "existing_t183_residual_datum_pair12_orbit_direction_selection_frontier"
ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
CURRENT_WORK_CONTRACT = "build_one_narrow_probe_for_a_genuinely_new_inversion_sensitive_source_side_provider_class_against_the_existing_f960_bridge_target"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1115, IN_F972, IN_P708, IN_P1060, IN_P1061, IN_F960]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1116",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1115 = load_json(IN_P1115)
    f972 = load_json(IN_F972)
    p708 = load_json(IN_P708)
    p1060 = load_json(IN_P1060)
    p1061 = load_json(IN_P1061)
    f960 = load_json(IN_F960)

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

    existing_f972_selector_source_subroute_confirmed = (
        p1115.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F966_NON_SAME_LANE_ROUTE_CONTRACT_TO_EXPORTED_STRICT_SELECTOR_SOURCE_FRONTIER_SUBROUTE_DECISION_AUDITED"
        and p1115.get("strongest_existing_exported_subroute") == PRIMARY_SUBROUTE
        and p1115.get("primary_existing_subroute_target") == PRIMARY_TARGET
        and f972.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F966_NON_SAME_LANE_ROUTE_CONTRACT_TO_EXPORTED_STRICT_SELECTOR_SOURCE_FRONTIER_SUBROUTE_PACKET_EXPORTED"
        and f972.get("packet_exported_on_current_repo_state") is True
        and f972.get("strongest_existing_exported_subroute") == PRIMARY_SUBROUTE
        and f972.get("primary_existing_subroute_target") == PRIMARY_TARGET
    )

    t173_residual_pair12_split_localized = (
        p708.get("status") == "PASS_T173_FRONTIER_DASHBOARD_READY"
        and p708.get("recommended_next_strict_target") == "T173"
        and p708.get("residual_datum_pair12_split_localized_as_opposite_orbit_directions") is True
        and p708.get("t183_residual_datum_pair12_orbit_direction_selection_bridge_exported") is False
        and p708.get("surviving_pair12_residual_datum_carrier_remains_selector_neutral") is True
    )

    existing_t183_frontier_confirmed = (
        p1060.get("status") == "PASS_CURRENT_STRICT_T173_T176_POST_STOP_SOURCE_SIDE_INPUT_LEG_LANE_TO_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_FRONTIER_ROUTE_DECISION_AUDITED"
        and p1060.get("existing_t183_frontier_named_and_open") is True
        and p1060.get("primary_continuation_route") == ACTIVE_FRONTIER
    )

    no_current_lawful_supplier_confirmed = (
        p1061.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_EXISTING_EXPORT_OR_NEAR_MISS_CANDIDATES_AUDITED_NO_LAWFUL_SUPPLIER_FOUND"
        and p1061.get("t183_target_name") == ACTIVE_BRIDGE
        and p1061.get("current_repo_already_exports_actual_t183_supplier") is False
    )

    exact_active_bridge_target_exported = (
        f960.get("status") == "F960_EXECUTED_CURRENT_STRICT_T173_T176_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
        and f960.get("target_object_id") == ACTIVE_BRIDGE
    )

    strongest_exact_active_frontier_is_existing_f960_t183_bridge = (
        existing_f972_selector_source_subroute_confirmed
        and t173_residual_pair12_split_localized
        and existing_t183_frontier_confirmed
        and no_current_lawful_supplier_confirmed
        and exact_active_bridge_target_exported
    )

    add_check("existing_f972_selector_source_subroute_confirmed", existing_f972_selector_source_subroute_confirmed, True, "P1115/F972 already freeze the strongest exported non-same-lane subroute as the strict selector-source frontier rejoining T173.")
    add_check("t173_residual_pair12_split_localized", t173_residual_pair12_split_localized, True, "P708 already localizes the active unresolved residual datum pair12 orbit-direction split under T173.")
    add_check("existing_t183_frontier_confirmed", existing_t183_frontier_confirmed, True, "P1060 already freezes the post-stop primary continuation route as the existing T183 frontier.")
    add_check("no_current_lawful_supplier_confirmed", no_current_lawful_supplier_confirmed, True, "P1061 already confirms that no lawful supplier is currently exported there.")
    add_check("exact_active_bridge_target_exported", exact_active_bridge_target_exported, True, "F960 already exports the exact active missing bridge target.")
    add_check("strongest_exact_active_frontier_is_existing_f960_t183_bridge", strongest_exact_active_frontier_is_existing_f960_t183_bridge, True, "Therefore the strongest exact active frontier beneath the live F972 selector-source subroute is the existing F960/T183 bridge frontier.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F972_SELECTOR_SOURCE_SUBROUTE_TO_EXISTING_F960_ACTIVE_BRIDGE_FRONTIER_ROUTE_DECISION_AUDITED"
        if not blocking and strongest_exact_active_frontier_is_existing_f960_t183_bridge
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F972_SELECTOR_SOURCE_SUBROUTE_TO_EXISTING_F960_ACTIVE_BRIDGE_FRONTIER_ROUTE_DECISION_AUDIT"
    )

    artifact = {
        "stage": "P1116",
        "status": status,
        "as_of": AS_OF,
        "strongest_exact_active_frontier": ACTIVE_FRONTIER,
        "strongest_exact_active_frontier_target": ACTIVE_BRIDGE,
        "current_primary_work_contract": CURRENT_WORK_CONTRACT,
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "strongest_exact_active_frontier": artifact["strongest_exact_active_frontier"],
        "strongest_exact_active_frontier_target": artifact["strongest_exact_active_frontier_target"],
        "current_primary_work_contract": artifact["current_primary_work_contract"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
