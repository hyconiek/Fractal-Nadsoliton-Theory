#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path


AS_OF = "2026-03-28"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1089 = GENERATED / "p1089_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_and_stop_audit_probe_summary.json"

OUT_JSON = GENERATED / "f964_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_stop_packet.json"
OUT_SUMMARY = GENERATED / "f964_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_stop_packet_summary.json"

P1089_STATUS = "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_STAGNATION_AND_STOP_AUDITED"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_P1089.exists():
        artifact = {
            "stage": "F964",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [rel(IN_P1089)],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1089 = load_json(IN_P1089)

    packet_exported_on_current_repo_state = (
        p1089.get("status") == P1089_STATUS
        and p1089.get("same_lane_stagnation_boundary_reached") is True
        and p1089.get("stop_condition_triggered") is True
    )

    status = (
        "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_STAGNATION_STOP_PACKET_EXPORTED"
        if packet_exported_on_current_repo_state
        else "FAIL_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_STAGNATION_STOP_PACKET"
    )

    packet = {
        "stage": "F964",
        "status": status,
        "as_of": AS_OF,
        "packet_name": "Xi_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_stop_packet_v1",
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "feeder_sharper_same_lane_stagnation_boundary_reached": p1089.get("same_lane_stagnation_boundary_reached"),
        "same_lane_deeper_witness_descent_disallowed_as_primary_move": p1089.get(
            "further_same_lane_sharper_witness_descent_is_not_honest_primary_move"
        ),
        "same_lane_stop_threshold_attempt_count": p1089.get("same_lane_stop_threshold_attempt_count"),
        "same_lane_exact_attempt_count": p1089.get("same_lane_exact_attempt_count"),
        "same_lane_sharper_target_count": p1089.get("same_lane_sharper_target_count"),
        "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route": p1089.get(
            "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"
        ),
        "current_primary_work_contract": "stop_same_lane_feeder_sharper_witness_descent_not_fake_progress",
        "no_false_pass": True,
    }
    write_json(OUT_JSON, packet)

    summary = {
        "stage": packet["stage"],
        "status": packet["status"],
        "as_of": packet["as_of"],
        "packet_exported_on_current_repo_state": packet["packet_exported_on_current_repo_state"],
        "same_lane_stagnation_boundary_reached": packet["feeder_sharper_same_lane_stagnation_boundary_reached"],
        "same_lane_deeper_witness_descent_disallowed_as_primary_move": packet[
            "same_lane_deeper_witness_descent_disallowed_as_primary_move"
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
