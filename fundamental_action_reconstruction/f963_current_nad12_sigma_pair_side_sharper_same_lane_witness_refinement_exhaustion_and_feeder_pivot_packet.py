#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1064 = GENERATED / "p1064_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_and_feeder_pivot_audit_probe_summary.json"
OUT_JSON = GENERATED / "f963_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_and_feeder_pivot_packet.json"
OUT_SUMMARY = GENERATED / "f963_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_and_feeder_pivot_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    if not IN_P1064.exists():
        artifact = {
            "stage": "F963",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P1064.relative_to(REPO))],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1064 = load_json(IN_P1064)
    admitted = bool(p1064.get("next_honest_move_is_noncyclic_pivot_to_feeder_support_side_witness_family"))

    status = (
        "F963_EXECUTED_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_EXHAUSTION_AND_FEEDER_PIVOT_PACKET_NO_FALSE_PASS"
        if admitted
        else "F963_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ROUTE_DECISION_STATE"
    )

    artifact = {
        "stage": "F963",
        "status": status,
        "as_of": AS_OF,
        "route_object_id": "Nad12SigmaPairSideSharperSameLaneExhaustionFeederPivot_v1",
        "same_lane_exhaustion_boundary_reached": bool(p1064.get("same_lane_exhaustion_boundary_reached")),
        "preferred_noncyclic_pivot_family": p1064.get("preferred_noncyclic_pivot_family"),
        "preferred_first_pivot_branch": p1064.get("preferred_first_pivot_branch"),
        "hard_limits": [
            "Does not export actual feeder-support-side provider support.",
            "Does not export actual pair-realization-side provider support.",
            "Does not export actual feeder support, theta export, pair population, residual bridge support, or loop break.",
            "Does not discharge QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "route_object_id": artifact["route_object_id"],
        "same_lane_exhaustion_boundary_reached": artifact["same_lane_exhaustion_boundary_reached"],
        "preferred_noncyclic_pivot_family": artifact["preferred_noncyclic_pivot_family"],
        "preferred_first_pivot_branch": artifact["preferred_first_pivot_branch"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
