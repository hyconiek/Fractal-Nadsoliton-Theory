#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1059 = GENERATED / "p1059_current_strict_t173_t176_source_side_input_leg_same_lane_stagnation_and_stop_audit_probe_summary.json"

OUT_JSON = GENERATED / "f958_current_strict_t173_t176_source_side_input_leg_same_lane_stagnation_stop_packet.json"
OUT_SUMMARY = GENERATED / "f958_current_strict_t173_t176_source_side_input_leg_same_lane_stagnation_stop_packet_summary.json"

P1059_STATUS = "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_SAME_LANE_STAGNATION_AND_STOP_AUDITED"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_P1059.exists():
        artifact = {
            "stage": "F958",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [rel(IN_P1059)],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1059 = load_json(IN_P1059)

    packet_exported_on_current_repo_state = (
        p1059.get("status") == P1059_STATUS
        and p1059.get("same_lane_stagnation_boundary_reached") is True
        and p1059.get("stop_condition_triggered") is True
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_SAME_LANE_STAGNATION_STOP_PACKET_EXPORTED"
        if packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_SAME_LANE_STAGNATION_STOP_PACKET"
    )

    packet = {
        "stage": "F958",
        "status": status,
        "as_of": AS_OF,
        "packet_name": "Xi_strict_t173_t176_source_side_input_leg_same_lane_stagnation_stop_packet_v1",
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "source_side_input_leg_same_lane_stagnation_boundary_reached": p1059.get(
            "same_lane_stagnation_boundary_reached"
        ),
        "same_lane_deeper_boundary_descent_disallowed_as_primary_move": p1059.get(
            "further_same_lane_descent_is_not_honest_primary_move"
        ),
        "same_lane_stop_threshold_attempt_count": p1059.get("same_lane_stop_threshold_attempt_count"),
        "same_lane_exact_attempt_count": p1059.get("same_lane_exact_attempt_count"),
        "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route": p1059.get(
            "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"
        ),
        "current_primary_work_contract": "stop_same_lane_source_side_input_leg_descent_not_fake_progress",
        "no_false_pass": True,
    }
    write_json(OUT_JSON, packet)

    summary = {
        "stage": packet["stage"],
        "status": packet["status"],
        "as_of": packet["as_of"],
        "packet_exported_on_current_repo_state": packet["packet_exported_on_current_repo_state"],
        "same_lane_stagnation_boundary_reached": packet["source_side_input_leg_same_lane_stagnation_boundary_reached"],
        "same_lane_deeper_boundary_descent_disallowed_as_primary_move": packet[
            "same_lane_deeper_boundary_descent_disallowed_as_primary_move"
        ],
        "same_lane_stop_threshold_attempt_count": packet["same_lane_stop_threshold_attempt_count"],
        "same_lane_exact_attempt_count": packet["same_lane_exact_attempt_count"],
        "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route": packet[
            "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
