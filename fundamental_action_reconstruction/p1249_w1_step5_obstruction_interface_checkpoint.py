#!/usr/bin/env python3
"""P1249: interface checkpoint between strict chain (C1.3) and QW-2191 obstruction."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1248", type=Path, default=GEN / "p1248_w1_step4_formal_chain_extension_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1249_w1_step5_obstruction_interface_summary.json")
    args = parser.parse_args()

    p1248 = json.loads(args.p1248.read_text(encoding="utf-8"))
    c1_3_supported = bool(p1248.get("c1_3_supported", False))

    coverage = {
        "strict_lane_internal_chain_consistency": c1_3_supported,
        "non_strict_dependency_excluded": c1_3_supported,
        "selector_uniqueness_obstruction_fully_discharged": False,
    }

    remaining_gaps = [
        "export a strict-core selector uniqueness witness that resolves QW-2191 without non-strict axiom route",
        "promote from local chain consistency to obstruction-level uniqueness discharge packet",
    ]

    out = {
        "packet": "P1249",
        "as_of": "2026-05-11",
        "source_step": "P1248_C1.3",
        "qw_2191_interface_coverage": coverage,
        "remaining_obstruction_gaps": remaining_gaps,
        "qw_2191_status": "OPEN",
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Obstruction interface checkpoint: chain progress acknowledged, QW-2191 not yet discharged.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1249] c1_3_supported={c1_3_supported} wrote {args.out}")


if __name__ == "__main__":
    main()
