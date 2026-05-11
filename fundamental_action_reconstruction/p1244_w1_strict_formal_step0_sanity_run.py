#!/usr/bin/env python3
"""P1244: first strict formal step-0 sanity run packet."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1243", type=Path, default=GEN / "p1243_w1_strict_statement_consistency_summary.json")
    parser.add_argument("--scaffold", type=Path, default=GEN / "p1239_w1_strict_formal_step0_scaffold.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1244_w1_strict_formal_step0_sanity_summary.json")
    args = parser.parse_args()

    p1243 = json.loads(args.p1243.read_text(encoding="utf-8"))
    scaffold = json.loads(args.scaffold.read_text(encoding="utf-8"))

    entry_ok = bool(p1243.get("strict_statement_consistent_with_trace", False))

    formal_packet = {
        "packet": "P1244_FORMAL_STEP0",
        "as_of": "2026-05-11",
        "formal_statement_v1": scaffold.get("w1_strict_formal_statement"),
        "assumptions_v1": scaffold.get("w1_strict_formal_assumptions", []),
        "derived_claims_v1": [
            "C1: strict-path discharge is locally admissible on strict lane under listed assumptions",
            "C2: non-strict axiom route is excluded from this strict-step0 derivation",
            "C3: theory closure status remains OPEN at this stage",
        ],
        "proof_stage": "STEP0_SANITY",
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
    }

    out = {
        "packet": "P1244",
        "as_of": "2026-05-11",
        "entry_condition_from_p1243": entry_ok,
        "formal_packet_emitted": entry_ok,
        "formal_packet": formal_packet if entry_ok else None,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "First strict formal sanity packet emitted from traced, consistent step-0 inputs.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1244] emitted={entry_ok} wrote {args.out}")


if __name__ == "__main__":
    main()
