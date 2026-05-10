#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1113 = GENERATED / "p1113_current_strict_t173_t176_minimal_onrd_ssl_rcw_stagnation_audit_probe_summary.json"

OUT_JSON = GENERATED / "f970_current_strict_t173_t176_minimal_onrd_ssl_rcw_stop_packet.json"
OUT_SUMMARY = GENERATED / "f970_current_strict_t173_t176_minimal_onrd_ssl_rcw_stop_packet_summary.json"

P1113_STATUS = "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STAGNATION_AND_STOP_AUDITED"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_P1113.exists():
        artifact = {
            "stage": "F970",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [rel(IN_P1113)],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1113 = load_json(IN_P1113)

    packet_exported_on_current_repo_state = (
        p1113.get("status") == P1113_STATUS
        and p1113.get("same_lane_stagnation_boundary_reached") is True
        and p1113.get("stop_condition_triggered") is True
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STAGNATION_STOP_PACKET_EXPORTED"
        if packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STAGNATION_STOP_PACKET"
    )

    packet = {
        "stage": "F970",
        "status": status,
        "as_of": AS_OF,
        "packet_name": "Xi_current_strict_t173_t176_minimal_onrd_sharper_same_lane_route_coherence_witness_refinement_stagnation_stop_packet_v1",
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "minimal_onrd_sharper_same_lane_stagnation_boundary_reached": p1113.get("same_lane_stagnation_boundary_reached"),
        "same_lane_deeper_route_coherence_witness_descent_disallowed_as_primary_move": p1113.get(
            "further_same_lane_sharper_route_coherence_witness_descent_is_not_honest_primary_move"
        ),
        "same_lane_stop_threshold_attempt_count": p1113.get("same_lane_stop_threshold_attempt_count"),
        "same_lane_exact_attempt_count": p1113.get("same_lane_exact_attempt_count"),
        "same_lane_sharper_target_count": p1113.get("same_lane_sharper_target_count"),
        "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route": p1113.get(
            "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"
        ),
        "current_primary_work_contract": "stop_same_lane_minimal_onrd_sharper_route_coherence_witness_descent_not_fake_progress",
        "no_false_pass": True,
    }
    write_json(OUT_JSON, packet)

    summary = {
        "stage": packet["stage"],
        "status": packet["status"],
        "as_of": packet["as_of"],
        "packet_exported_on_current_repo_state": packet["packet_exported_on_current_repo_state"],
        "same_lane_stagnation_boundary_reached": packet["minimal_onrd_sharper_same_lane_stagnation_boundary_reached"],
        "same_lane_deeper_route_coherence_witness_descent_disallowed_as_primary_move": packet[
            "same_lane_deeper_route_coherence_witness_descent_disallowed_as_primary_move"
        ],
        "same_lane_stop_threshold_attempt_count": packet["same_lane_stop_threshold_attempt_count"],
        "same_lane_exact_attempt_count": packet["same_lane_exact_attempt_count"],
        "same_lane_sharper_target_count": packet["same_lane_sharper_target_count"],
        "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route": packet[
            "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
