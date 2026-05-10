#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1116 = GENERATED / "p1116_current_strict_t173_t176_f972_to_f960_active_bridge_route_audit_probe_summary.json"
IN_F973 = GENERATED / "f973_current_strict_t173_t176_f972_to_f960_active_bridge_route_packet_summary.json"
IN_P1063 = GENERATED / "p1063_current_strict_t173_t176_post_f961_negative_3_cycle_reference_to_existing_t216_t218_pair12_provider_frontier_route_decision_audit_probe_summary.json"
IN_F962 = GENERATED / "f962_current_strict_t173_t176_post_f961_negative_3_cycle_reference_to_existing_t216_t218_pair12_provider_frontier_route_decision_packet_summary.json"
IN_P982 = GENERATED / "p982_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_audit_probe_summary.json"
IN_F949 = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_summary.json"
IN_P1091 = GENERATED / "p1091_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_decision_audit_probe_summary.json"
IN_F966 = GENERATED / "f966_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_packet_summary.json"

OUT_JSON = GENERATED / "p1117_current_strict_t173_t176_f973_f960_no_exported_live_nonsamelane_provider_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1117_current_strict_t173_t176_f973_f960_no_exported_live_nonsamelane_provider_audit_probe_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
ACTIVE_FRONTIER = "existing_t183_residual_datum_pair12_orbit_direction_selection_frontier"
NEW_WORK_CONTRACT = "build_one_genuinely_new_narrow_probe_against_existing_f960_bridge_target"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1116, IN_F973, IN_P1063, IN_F962, IN_P982, IN_F949, IN_P1091, IN_F966]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1117",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1116 = load_json(IN_P1116)
    f973 = load_json(IN_F973)
    p1063 = load_json(IN_P1063)
    f962 = load_json(IN_F962)
    p982 = load_json(IN_P982)
    f949 = load_json(IN_F949)
    p1091 = load_json(IN_P1091)
    f966 = load_json(IN_F966)

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

    active_f960_bridge_frontier_confirmed = (
        p1116.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F972_SELECTOR_SOURCE_SUBROUTE_TO_EXISTING_F960_ACTIVE_BRIDGE_FRONTIER_ROUTE_DECISION_AUDITED"
        and p1116.get("strongest_exact_active_frontier") == ACTIVE_FRONTIER
        and p1116.get("strongest_exact_active_frontier_target") == ACTIVE_BRIDGE
        and f973.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F972_SELECTOR_SOURCE_SUBROUTE_TO_EXISTING_F960_ACTIVE_BRIDGE_FRONTIER_ROUTE_PACKET_EXPORTED"
        and f973.get("packet_exported_on_current_repo_state") is True
        and f973.get("strongest_exact_active_frontier") == ACTIVE_FRONTIER
        and f973.get("strongest_exact_active_frontier_target") == ACTIVE_BRIDGE
    )

    old_pair12_provider_frontier_exists_but_is_nonprimary = (
        p1063.get("status") == "PASS_CURRENT_STRICT_T173_T176_POST_F961_NEGATIVE_3_CYCLE_REFERENCE_TO_EXISTING_T216_T218_PAIR12_PROVIDER_FRONTIER_ROUTE_DECISION_AUDITED"
        and p1063.get("rejoin_to_existing_t216_t218_frontier") is True
        and p1063.get("primary_continuation_target") == "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_v1"
        and p1063.get("exact_interface_frontier_target") == "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_v1"
        and f962.get("status") == "F962_EXECUTED_CURRENT_STRICT_T173_T176_POST_F961_NEGATIVE_3_CYCLE_REFERENCE_TO_EXISTING_T216_T218_PAIR12_PROVIDER_FRONTIER_ROUTE_DECISION_PACKET_NO_FALSE_PASS"
        and f962.get("rejoin_to_existing_t216_t218_frontier") is True
    )

    pair12_same_lane_reentry_disallowed = (
        p982.get("status") == "PASS_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_AUDITED"
        and p982.get("same_lane_exhaustion_boundary_reached") is True
        and p982.get("further_same_lane_t274_style_descent_is_not_honest_primary_move") is True
        and f949.get("status") == "PASS_CURRENT_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_PACKET_EXPORTED"
        and f949.get("same_lane_exhaustion_boundary_reached") is True
    )

    live_contract_requires_genuinely_new_non_same_lane_or_new_provider_class = (
        p1091.get("status") == "PASS_CURRENT_STRICT_T173_T176_POST_F965_FAILURE_MAP_TO_EXPORTED_NONCYCLIC_PROVIDER_SPLIT_NON_SAME_LANE_UPGRADE_ROUTE_DECISION_AUDITED"
        and p1091.get("pair12_entry_same_lane_reentry_disallowed_as_primary_move") is True
        and p1091.get("active_missing_bridge") == ACTIVE_BRIDGE
        and f966.get("status") == "F966_EXECUTED_CURRENT_STRICT_T173_T176_POST_F965_FAILURE_MAP_TO_EXPORTED_NONCYCLIC_PROVIDER_SPLIT_NON_SAME_LANE_UPGRADE_ROUTE_PACKET_NO_FALSE_PASS"
        and f966.get("allowed_next_move_contract") == "search_one_genuinely_new_non_same_lane_upgrade_route_within_exported_noncyclic_provider_split_family_or_one_genuinely_new_inversion_sensitive_source_side_provider_class"
    )

    already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present = False

    add_check("active_f960_bridge_frontier_confirmed", active_f960_bridge_frontier_confirmed, True, "P1116/F973 already freeze the active exact bridge frontier as F960/T183.")
    add_check("old_pair12_provider_frontier_exists_but_is_nonprimary", old_pair12_provider_frontier_exists_but_is_nonprimary, True, "P1063/F962 already show one old exact pair12 provider frontier inside T216/T218.")
    add_check("pair12_same_lane_reentry_disallowed", pair12_same_lane_reentry_disallowed, True, "P982/F949 already freeze that reopening the old pair12 entry same-lane descent is not an honest primary move.")
    add_check("live_contract_requires_genuinely_new_non_same_lane_or_new_provider_class", live_contract_requires_genuinely_new_non_same_lane_or_new_provider_class, True, "P1091/F966 already freeze that the admissible next-move contract is one genuinely new non-same-lane upgrade route or one genuinely new inversion-sensitive source-side provider class.")
    add_check("already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present", already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present, False, "Therefore the current repo still exports no already-live non-same-lane inversion-sensitive source-side provider candidate beneath F960.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F973_F960_NO_ALREADY_EXPORTED_LIVE_NON_SAME_LANE_INVERSION_SENSITIVE_SOURCE_SIDE_PROVIDER_CANDIDATE_AUDITED"
        if not blocking else
        "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F973_F960_NO_ALREADY_EXPORTED_LIVE_NON_SAME_LANE_INVERSION_SENSITIVE_SOURCE_SIDE_PROVIDER_CANDIDATE_AUDIT"
    )

    artifact = {
        "stage": "P1117",
        "status": status,
        "as_of": AS_OF,
        "already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present": already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "current_primary_work_contract": NEW_WORK_CONTRACT,
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present": artifact["already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "current_primary_work_contract": artifact["current_primary_work_contract"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
