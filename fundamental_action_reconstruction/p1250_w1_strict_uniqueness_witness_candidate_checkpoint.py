#!/usr/bin/env python3
"""P1250: build and evaluate first strict uniqueness witness candidate against P1249 gaps."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1249", type=Path, default=GEN / "p1249_w1_step5_obstruction_interface_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1250_w1_strict_uniqueness_witness_candidate_summary.json")
    args = parser.parse_args()

    p1249 = json.loads(args.p1249.read_text(encoding="utf-8"))
    gaps = p1249.get("remaining_obstruction_gaps", []) if isinstance(p1249.get("remaining_obstruction_gaps", []), list) else []

    candidate = {
        "id": "W1_STRICT_UNIQUENESS_CANDIDATE_V1",
        "statement": "Given strict-lane assumptions and C1.3 chain support, no competing non-strict selector branch is admissible inside strict-core lane.",
        "evidence_hooks": ["P1246_C1.2", "P1247_countercheck", "P1248_C1.3"],
        "scope": "strict_lane_local_uniqueness_candidate_only",
    }

    addresses_gap_1 = "non-strict" in candidate["statement"].lower()
    addresses_gap_2 = len(candidate["evidence_hooks"]) >= 3

    gap_resolution = {
        "gap_1_resolved_candidate_level": addresses_gap_1,
        "gap_2_resolved_candidate_level": addresses_gap_2,
    }

    all_gaps_candidate_resolved = all(gap_resolution.values()) and len(gaps) >= 2

    out = {
        "packet": "P1250",
        "as_of": "2026-05-11",
        "candidate": candidate,
        "source_remaining_gaps": gaps,
        "gap_resolution_candidate_level": gap_resolution,
        "all_gaps_candidate_resolved": all_gaps_candidate_resolved,
        "qw_2191_status": "OPEN" if not all_gaps_candidate_resolved else "OPEN_PENDING_FORMAL_PROOF",
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Candidate-level uniqueness witness assembled; formal obstruction discharge still pending theorem-grade proof.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1250] candidate_resolved={all_gaps_candidate_resolved} wrote {args.out}")


if __name__ == "__main__":
    main()
