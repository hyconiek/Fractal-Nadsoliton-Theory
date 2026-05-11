#!/usr/bin/env python3
"""P1233: strict-lane readiness gate from comparative W1 checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1232", type=Path, default=GEN / "p1232_w1_strict_vs_nonstrict_comparative_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1233_w1_strict_lane_readiness_gate_summary.json")
    args = parser.parse_args()

    p1232 = json.loads(args.p1232.read_text(encoding="utf-8"))
    rows = p1232.get("comparison", []) if isinstance(p1232.get("comparison", []), list) else []

    strict_row = next((r for r in rows if isinstance(r, dict) and r.get("lane") == "STRICT_CANDIDATE"), {})

    cond1 = strict_row.get("discharge_mode") == "STRICT_PATH_DISCHARGE"
    cond2 = bool(strict_row.get("strict_core_closure_eligibility", False))
    cond3 = strict_row.get("symmetry_state") == "UNBROKEN_IN_STRICT_CORE"

    gate_pass = cond1 and cond2 and cond3

    out = {
        "packet": "P1233",
        "as_of": "2026-05-11",
        "strict_lane_conditions": {
            "discharge_mode_is_strict_path": cond1,
            "strict_core_closure_eligibility_true": cond2,
            "symmetry_state_unbroken_in_strict_core": cond3,
        },
        "strict_lane_readiness_gate_pass": gate_pass,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Readiness gate only; pass does not imply global closure.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1233] gate_pass={gate_pass} wrote {args.out}")


if __name__ == "__main__":
    main()
