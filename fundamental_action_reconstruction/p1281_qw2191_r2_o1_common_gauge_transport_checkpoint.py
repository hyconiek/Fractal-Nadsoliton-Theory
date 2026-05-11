#!/usr/bin/env python3
"""P1281: R2.O1 common-gauge transport checkpoint for QW-2191 bounds."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1280", type=Path, default=GEN / "p1280_qw2191_r2_bound_transport_checkpoint_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1281_qw2191_r2_o1_common_gauge_transport_summary.json")
    args = parser.parse_args()

    p1280 = _read(args.p1280)
    if p1280.get("next_priority") != "R2_O1_COMMON_GAUGE_TRANSPORT":
        raise SystemExit("P1281 requires next_priority=R2_O1_COMMON_GAUGE_TRANSPORT from P1280.")
    if p1280.get("closure_policy", {}).get("global_qw2191_closure_allowed", True):
        raise SystemExit("P1281 requires unresolved global QW-2191 closure policy from P1280.")

    gauges = [
        "G_local_sector_A",
        "G_local_sector_B",
        "G_local_sector_C",
    ]
    target_gauge = "G_qw2191_common_reference_v1"
    transport_maps = [
        {"from": g, "to": target_gauge, "map_id": f"tau_{idx+1}", "status": "DECLARED"}
        for idx, g in enumerate(gauges)
    ]

    out = {
        "packet": "P1281",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1280": str(args.p1280)},
        "r2_o1": {
            "target_gauge": target_gauge,
            "source_gauges": gauges,
            "transport_maps": transport_maps,
            "consistency_certificate": "PENDING",
            "status": "PARTIAL_DISCHARGE",
        },
        "residual_risk": {
            "mismatch_control_lemma_required": True,
            "certificate_required": True,
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
            "bridge_decision_required": True,
        },
        "next_priority": "R2_O2_MISMATCH_CONTROL_LEMMA",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1281] wrote {args.out}; maps={len(transport_maps)} status=PARTIAL_DISCHARGE")


if __name__ == "__main__":
    main()
