#!/usr/bin/env python3
"""P1221: continue workflow after local W4 discharge checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1210", type=Path, default=GEN / "p1210_controlled_w4_symbolic_attempt_summary.json")
    parser.add_argument("--p1192", type=Path, default=GEN / "p1192_selector_premise_witness_map_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1221_post_w4_discharge_continuation_summary.json")
    args = parser.parse_args()

    p1210 = json.loads(args.p1210.read_text(encoding="utf-8"))
    p1192 = json.loads(args.p1192.read_text(encoding="utf-8"))

    w4_discharged = bool(p1210.get("w4_discharge_pass", False))
    open_count = int(p1192.get("open_count", 0)) if isinstance(p1192.get("open_count", 0), int) else 0

    out = {
        "packet": "P1221",
        "as_of": "2026-05-11",
        "w4_local_discharge_checkpoint": w4_discharged,
        "w4_status": "LOCALLY_DISCHARGED" if w4_discharged else "OPEN",
        "remaining_open_witness_count_reported": max(0, open_count - (1 if w4_discharged else 0)),
        "next_target": "ATTACK_NEXT_OPEN_WITNESS_OBLIGATION" if w4_discharged else "KEEP_W4_ROUTE_ACTIVE",
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Post-W4 local continuation checkpoint; no global closure claim.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1221] w4_local_discharge_checkpoint={w4_discharged} wrote {args.out}")


if __name__ == "__main__":
    main()
