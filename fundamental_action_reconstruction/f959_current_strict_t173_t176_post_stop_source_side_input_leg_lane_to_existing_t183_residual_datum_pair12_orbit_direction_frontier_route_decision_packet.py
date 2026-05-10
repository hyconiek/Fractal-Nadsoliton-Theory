#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1060 = GENERATED / "p1060_current_strict_t173_t176_post_stop_source_side_input_leg_lane_to_existing_t183_residual_datum_pair12_orbit_direction_frontier_route_decision_audit_probe_summary.json"
IN_F959 = ROOT / "F959_CURRENT_STRICT_T173_T176_POST_STOP_SOURCE_SIDE_INPUT_LEG_LANE_TO_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_FRONTIER_ROUTE_DECISION_PACKET.md"

OUT_JSON = GENERATED / "f959_current_strict_t173_t176_post_stop_source_side_input_leg_lane_to_existing_t183_residual_datum_pair12_orbit_direction_frontier_route_decision_packet.json"
OUT_SUMMARY = GENERATED / "f959_current_strict_t173_t176_post_stop_source_side_input_leg_lane_to_existing_t183_residual_datum_pair12_orbit_direction_frontier_route_decision_packet_summary.json"


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

    prerequisites = [IN_P1060, IN_F959]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F959",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1060 = load_json(IN_P1060)
    f959_text = load_text(IN_F959)

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

    p1060_route_decision_passed = (
        p1060.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_POST_STOP_SOURCE_SIDE_INPUT_LEG_LANE_TO_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_FRONTIER_ROUTE_DECISION_AUDITED"
        and p1060.get("primary_continuation_route")
        == "existing_t183_residual_datum_pair12_orbit_direction_selection_frontier"
    )

    f959_packet_shape_frozen = all(
        needle in f959_text
        for needle in [
            "Xi_strict_t173_t176_post_stop_source_side_input_leg_lane_to_existing_t183_residual_datum_pair12_orbit_direction_frontier_route_decision_packet_v1",
            "stopped_source_side_input_leg_lane_confirmed := yes",
            "route_requires_non_same_lane_upgrade := yes",
            "existing_t183_frontier_named_and_open := yes",
            "direction_free_shannon_t184_negative_boundary_confirmed := yes",
            "t185_w_break_family_not_admitted_as_active_primary_move := yes",
            "primary_continuation_route := existing_t183_residual_datum_pair12_orbit_direction_selection_frontier",
            "current_primary_work_contract := leave_stopped_source_side_input_leg_same_lane_and_do_not_reactivate_negative_t184_or_t185_continuation_families",
        ]
    )

    packet_exported_on_current_repo_state = p1060_route_decision_passed and f959_packet_shape_frozen

    add_check(
        "p1060_route_decision_passed",
        p1060_route_decision_passed,
        True,
        "P1060 already freezes the post-stop route-decision audit positively.",
    )
    add_check(
        "f959_packet_shape_frozen",
        f959_packet_shape_frozen,
        True,
        "F959 freezes the route-decision packet shape explicitly.",
    )
    add_check(
        "packet_exported_on_current_repo_state",
        packet_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one honest post-stop route-decision packet to the existing T183 frontier.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_POST_STOP_SOURCE_SIDE_INPUT_LEG_LANE_TO_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_FRONTIER_ROUTE_DECISION_PACKET_EXPORTED"
        if not blocking and packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_POST_STOP_SOURCE_SIDE_INPUT_LEG_LANE_TO_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_FRONTIER_ROUTE_DECISION_PACKET"
    )

    artifact = {
        "stage": "F959",
        "status": status,
        "as_of": AS_OF,
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "primary_continuation_route": "existing_t183_residual_datum_pair12_orbit_direction_selection_frontier",
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "packet_exported_on_current_repo_state": artifact["packet_exported_on_current_repo_state"],
        "primary_continuation_route": artifact["primary_continuation_route"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
