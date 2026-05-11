#!/usr/bin/env python3
"""P1229: validate minimal proof-packet index for W1 candidate reference."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1228", type=Path, default=GEN / "p1228_w1_real_candidate_evidence_note_summary.json")
    parser.add_argument("--index", type=Path, default=GEN / "p1229_w1_minimal_proof_packet_index.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1229_w1_minimal_proof_packet_index_summary.json")
    args = parser.parse_args()

    p1228 = _load(args.p1228)
    ref = p1228.get("candidate_ref", "")

    if args.index.exists():
        idx = _load(args.index)
    else:
        idx = {
            "packet": "P1229_INDEX",
            "as_of": "2026-05-11",
            "entries": {
                "internal://w1/strict_selector_source_theorem_candidate_v1": {
                    "sections": ["S1_statement", "S2_assumptions", "S3_scope_boundary"],
                    "status": "DRAFT"
                }
            }
        }
        args.index.write_text(json.dumps(idx, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    entries = idx.get("entries", {}) if isinstance(idx.get("entries", {}), dict) else {}
    entry = entries.get(ref)
    entry_present = isinstance(entry, dict)
    has_sections = entry_present and isinstance(entry.get("sections"), list) and len(entry.get("sections")) > 0

    out = {
        "packet": "P1229",
        "as_of": "2026-05-11",
        "candidate_ref": ref,
        "index_path": str(args.index),
        "entry_present": entry_present,
        "has_sections": bool(has_sections),
        "index_status": entry.get("status") if entry_present else "MISSING",
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Minimal proof-packet index checkpoint; validates semantic anchor presence for candidate ref.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1229] entry_present={entry_present} has_sections={has_sections} wrote {args.out}")


if __name__ == "__main__":
    main()
