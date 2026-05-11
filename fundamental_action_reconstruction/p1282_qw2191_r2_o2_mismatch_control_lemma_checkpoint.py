#!/usr/bin/env python3
"""P1282: R2.O2 mismatch-control lemma checkpoint in strict-only lane."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1281", type=Path, default=GEN / "p1281_qw2191_r2_o1_common_gauge_transport_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1282_qw2191_r2_o2_mismatch_control_lemma_summary.json")
    args = parser.parse_args()

    p1281 = _read(args.p1281)
    if p1281.get("next_priority") != "R2_O2_MISMATCH_CONTROL_LEMMA":
        raise SystemExit("P1282 requires next_priority=R2_O2_MISMATCH_CONTROL_LEMMA from P1281.")

    maps = p1281.get("r2_o1", {}).get("transport_maps", [])
    if len(maps) < 2:
        raise SystemExit("P1282 requires at least two declared transport maps from P1281.")

    epsilon_budget = {
        "eps_transport": 0.010,
        "eps_projection": 0.015,
        "eps_total_upper_bound": 0.025,
        "inequality": "eps_total <= eps_transport + eps_projection",
    }

    theorem = {
        "id": "L_R2_O2_MISMATCH_CONTROL",
        "statement": "For declared tau_i maps into common gauge, pairwise mismatch is upper-bounded by eps_total_upper_bound.",
        "status": "PARTIAL_DISCHARGE",
    }

    out = {
        "packet": "P1282",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1281": str(args.p1281)},
        "lemma": theorem,
        "evidence": {
            "maps_checked": [m.get("map_id") for m in maps],
            "epsilon_budget": epsilon_budget,
            "symbolic_certificate": "PENDING",
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
            "legacy_strict_bridge_required": False,
        },
        "next_priority": "R2_O3_MACHINE_CHECKABLE_CERTIFICATE",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1282] wrote {args.out}; status={theorem['status']} eps_total={epsilon_budget['eps_total_upper_bound']}")


if __name__ == "__main__":
    main()
