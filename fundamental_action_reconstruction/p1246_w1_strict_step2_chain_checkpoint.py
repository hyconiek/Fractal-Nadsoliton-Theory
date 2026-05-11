#!/usr/bin/env python3
"""P1246: strict Step-2 chain checkpoint deriving C1.2 from C1.1."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1245", type=Path, default=GEN / "p1245_w1_strict_step1_inference_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1246_w1_strict_step2_chain_summary.json")
    args = parser.parse_args()

    p1245 = json.loads(args.p1245.read_text(encoding="utf-8"))
    c1_1 = p1245.get("derived_inference", {}) if isinstance(p1245.get("derived_inference"), dict) else {}

    c1_1_supported = bool(c1_1.get("supported", False))
    inputs = p1245.get("inference_inputs", {}) if isinstance(p1245.get("inference_inputs"), dict) else {}
    has_ref = bool(inputs.get("has_ref_assumption", False))
    has_sym = bool(inputs.get("has_unbroken_assumption", False))
    strict_lane_stmt = bool(inputs.get("statement_is_strict_lane", False))

    c1_2_supported = c1_1_supported and has_ref and has_sym and strict_lane_stmt
    c1_2_text = (
        "C1.2: Given C1.1 and preserved strict-lane premises, the strict witness lane remains chain-consistent "
        "through Step-2 without requiring non-strict symmetry-breaking input."
    ) if c1_2_supported else "C1.2 not derivable from current Step-1 state"

    out = {
        "packet": "P1246",
        "as_of": "2026-05-11",
        "chain_rule": "C1.1 + (A_ref & A_sym & S_strict_lane) -> C1.2",
        "c1_1_supported": c1_1_supported,
        "c1_2_supported": c1_2_supported,
        "derived_inference": {
            "id": "C1.2",
            "text": c1_2_text,
        },
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Strict Step-2 chain checkpoint from Step-1 inference.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1246] C1.2_supported={c1_2_supported} wrote {args.out}")


if __name__ == "__main__":
    main()
