#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F958 = GENERATED / "f958_current_strict_t173_t176_source_side_input_leg_same_lane_stagnation_stop_packet_summary.json"
IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"
IN_P728 = GENERATED / "p728_current_strict_t182_residual_datum_source_side_boundary_shielded_sublane_reduction_audit_probe_summary.json"
IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P730 = GENERATED / "p730_current_strict_t184_direction_free_shannon_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_N725 = ROOT / "N725_CURRENT_STRICT_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_NONEXPORT_BOUNDARY_THEOREM.md"

OUT_JSON = GENERATED / "p1060_current_strict_t173_t176_post_stop_source_side_input_leg_lane_to_existing_t183_residual_datum_pair12_orbit_direction_frontier_route_decision_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1060_current_strict_t173_t176_post_stop_source_side_input_leg_lane_to_existing_t183_residual_datum_pair12_orbit_direction_frontier_route_decision_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F958, IN_P708, IN_P728, IN_P729, IN_P730, IN_P731, IN_N725]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1060",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f958 = load_json(IN_F958)
    p708 = load_json(IN_P708)
    p728 = load_json(IN_P728)
    p729 = load_json(IN_P729)
    p730 = load_json(IN_P730)
    p731 = load_json(IN_P731)
    n725_text = load_text(IN_N725)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    stopped_source_side_input_leg_lane_confirmed = (
        f958.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_SAME_LANE_STAGNATION_STOP_PACKET_EXPORTED"
        and f958.get("packet_exported_on_current_repo_state") is True
        and f958.get("same_lane_stagnation_boundary_reached") is True
        and f958.get("same_lane_deeper_boundary_descent_disallowed_as_primary_move") is True
        and f958.get("restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route") is True
    )

    route_requires_non_same_lane_upgrade = (
        stopped_source_side_input_leg_lane_confirmed
        and f958.get("same_lane_exact_attempt_count") == 3
    )

    existing_t183_frontier_named_and_open = (
        p728.get("status")
        == "PARTIAL_RESIDUAL_DATUM_SOURCE_SIDE_REDUCES_POSITIVE_CORRIDOR_TO_BOUNDARY_SHIELDED_SUBLANE_ONLY"
        and p728.get("residual_datum_source_side_supported_positive_charts") == ["pair1", "pair2"]
        and p728.get("current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane")
        is True
        and p729.get("status")
        == "PASS_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p729.get("t183_target_name")
        == "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
        and p729.get("t183_target_exported_on_current_repo_state") is False
        and p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions") is True
        and p708.get("t183_residual_datum_pair12_orbit_direction_selection_bridge_exported") is False
        and "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1" in n725_text
        and "the missing step is now finer" in n725_text.lower()
    )

    direction_free_shannon_t184_negative_boundary_confirmed = (
        p730.get("status")
        == "PASS_DIRECTION_FREE_SHANNON_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p730.get("t184_target_exported_on_current_repo_state") is False
        and p730.get("current_direction_free_shannon_lane_already_exports_pair1_pair2_o2_to_z2_cuts")
        is True
        and p730.get("current_direction_free_shannon_lane_selects_pair12_orbit_direction_branch") is False
        and p708.get("t184_direction_free_shannon_pair12_orbit_direction_selection_bridge_exported")
        is False
    )

    t185_w_break_family_not_admitted_as_active_primary_move = (
        p731.get("status")
        == "PASS_W_BREAK_WITNESS_PAYLOAD_PAIR12_ORBIT_DIRECTION_PROMOTION_BRIDGE_NONEXPORT_AUDITED"
        and p731.get("t185_target_exported_on_current_repo_state") is False
        and p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches") is True
        and p708.get("current_pair12_witness_split_current_exported_continuation_family_named_members_all_real")
        is True
        and p708.get("current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative")
        is True
        and p708.get("same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move")
        is True
    )

    primary_continuation_route_is_existing_t183_frontier = (
        route_requires_non_same_lane_upgrade
        and existing_t183_frontier_named_and_open
        and direction_free_shannon_t184_negative_boundary_confirmed
        and t185_w_break_family_not_admitted_as_active_primary_move
    )

    add_check(
        "stopped_source_side_input_leg_lane_confirmed",
        stopped_source_side_input_leg_lane_confirmed,
        True,
        "F958 already stops the source-side input-leg same lane as a primary strategy.",
    )
    add_check(
        "route_requires_non_same_lane_upgrade",
        route_requires_non_same_lane_upgrade,
        True,
        "The exported stop packet already requires a non-same-lane upgrade route.",
    )
    add_check(
        "existing_t183_frontier_named_and_open",
        existing_t183_frontier_named_and_open,
        True,
        "The sharper residual-datum pair1/pair2 orbit-direction frontier is already theorem-localized and still open.",
    )
    add_check(
        "direction_free_shannon_t184_negative_boundary_confirmed",
        direction_free_shannon_t184_negative_boundary_confirmed,
        True,
        "The current direction-free Shannon lane is already packaged as insufficient for that surviving pair1/pair2 split.",
    )
    add_check(
        "t185_w_break_family_not_admitted_as_active_primary_move",
        t185_w_break_family_not_admitted_as_active_primary_move,
        True,
        "The current witness-payload continuation family is explicitly no longer admitted as the active primary T173 move.",
    )
    add_check(
        "primary_continuation_route_is_existing_t183_frontier",
        primary_continuation_route_is_existing_t183_frontier,
        True,
        "Therefore the strongest honest post-stop continuation is to rejoin the existing T183 frontier.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_POST_STOP_SOURCE_SIDE_INPUT_LEG_LANE_TO_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_FRONTIER_ROUTE_DECISION_AUDITED"
        if not blocking and primary_continuation_route_is_existing_t183_frontier
        else "FAIL_CURRENT_STRICT_T173_T176_POST_STOP_SOURCE_SIDE_INPUT_LEG_LANE_TO_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_FRONTIER_ROUTE_DECISION_AUDIT"
    )

    artifact = {
        "stage": "P1060",
        "status": status,
        "as_of": AS_OF,
        "stopped_source_side_input_leg_lane_confirmed": stopped_source_side_input_leg_lane_confirmed,
        "route_requires_non_same_lane_upgrade": route_requires_non_same_lane_upgrade,
        "existing_t183_frontier_named_and_open": existing_t183_frontier_named_and_open,
        "direction_free_shannon_t184_negative_boundary_confirmed": direction_free_shannon_t184_negative_boundary_confirmed,
        "t185_w_break_family_not_admitted_as_active_primary_move": t185_w_break_family_not_admitted_as_active_primary_move,
        "primary_continuation_route": (
            "existing_t183_residual_datum_pair12_orbit_direction_selection_frontier"
            if primary_continuation_route_is_existing_t183_frontier
            else None
        ),
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "stopped_source_side_input_leg_lane_confirmed": artifact["stopped_source_side_input_leg_lane_confirmed"],
        "route_requires_non_same_lane_upgrade": artifact["route_requires_non_same_lane_upgrade"],
        "existing_t183_frontier_named_and_open": artifact["existing_t183_frontier_named_and_open"],
        "direction_free_shannon_t184_negative_boundary_confirmed": artifact[
            "direction_free_shannon_t184_negative_boundary_confirmed"
        ],
        "t185_w_break_family_not_admitted_as_active_primary_move": artifact[
            "t185_w_break_family_not_admitted_as_active_primary_move"
        ],
        "primary_continuation_route": artifact["primary_continuation_route"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
