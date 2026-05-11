#!/usr/bin/env python3
"""P1245: first strict Step-1 inference from P1244 formal step-0 packet."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1244", type=Path, default=GEN / "p1244_w1_strict_formal_step0_sanity_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1245_w1_strict_step1_inference_summary.json")
    args = parser.parse_args()

    p1244 = json.loads(args.p1244.read_text(encoding="utf-8"))
    formal_packet = p1244.get("formal_packet") if isinstance(p1244.get("formal_packet"), dict) else {}

    assumptions = formal_packet.get("assumptions_v1", []) if isinstance(formal_packet.get("assumptions_v1", []), list) else []
    statement = str(formal_packet.get("formal_statement_v1", ""))

    has_ref_assumption = any("strict_selector_source_theorem_exported=true" in a for a in assumptions)
    has_unbroken_assumption = any("UNBROKEN_IN_STRICT_CORE" in a for a in assumptions)
    has_no_global_assumption = any("theory_closure_status remains OPEN" in a for a in assumptions)
    statement_is_strict_lane = "strict lane" in statement.lower()

    c1_1_support = has_ref_assumption and has_unbroken_assumption and statement_is_strict_lane
    c1_1 = (
        "C1.1: Under strict-lane assumptions A_ref + A_sym and strict-lane statement S1, "
        "local W1 strict-path admissibility remains internally consistent at Step-1."
    ) if c1_1_support else "C1.1 not derivable from current step-0 packet"

    out = {
        "packet": "P1245",
        "as_of": "2026-05-11",
        "step0_packet_present": bool(formal_packet),
        "inference_inputs": {
            "has_ref_assumption": has_ref_assumption,
            "has_unbroken_assumption": has_unbroken_assumption,
            "has_no_global_assumption": has_no_global_assumption,
            "statement_is_strict_lane": statement_is_strict_lane,
        },
        "derived_inference": {
            "id": "C1.1",
            "supported": c1_1_support,
            "text": c1_1,
        },
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "First strict Step-1 inference checkpoint from P1244 formal packet.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1245] C1.1_supported={c1_1_support} wrote {args.out}")


if __name__ == "__main__":
    main()
