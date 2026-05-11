#!/usr/bin/env python3
"""P1228: attach minimal semantic evidence note to W1 real-candidate reference run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1227", type=Path, default=GEN / "p1227_w1_real_candidate_run_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1228_w1_real_candidate_evidence_note_summary.json")
    args = parser.parse_args()

    p1227 = json.loads(args.p1227.read_text(encoding="utf-8"))
    ref = p1227.get("candidate_ref", "")

    out = {
        "packet": "P1228",
        "as_of": "2026-05-11",
        "candidate_ref": ref,
        "evidence_scope_note": "reference marks a strict-path candidate handle for W1 continuity checks; it is not yet a full exported selector theorem proof packet",
        "assumption_boundary": "no global closure claim; QW-2191 remains formally open at theory level",
        "continuity_status": {
            "w1_witness_status": p1227.get("w1_witness_status"),
            "w1_discharge_mode": p1227.get("w1_discharge_mode"),
            "strict_binding_ok": p1227.get("strict_binding_ok"),
            "strict_reference_tier": p1227.get("strict_reference_tier"),
        },
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Semantic annotation checkpoint over P1227 result.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1228] wrote {args.out}")


if __name__ == "__main__":
    main()
