#!/usr/bin/env python3
"""P1234: freeze-note for strict-lane artifact set after readiness gate pass."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1233", type=Path, default=GEN / "p1233_w1_strict_lane_readiness_gate_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1234_w1_strict_lane_artifact_freeze_note_summary.json")
    args = parser.parse_args()

    p1233 = json.loads(args.p1233.read_text(encoding="utf-8"))
    gate_pass = bool(p1233.get("strict_lane_readiness_gate_pass", False))

    artifact_set = [
        "generated/p1227_w1_real_candidate_run_summary.json",
        "generated/p1228_w1_real_candidate_evidence_note_summary.json",
        "generated/p1229_w1_minimal_proof_packet_index.json",
        "generated/p1230_w1_section_completeness_and_symmetry_summary.json",
        "generated/p1232_w1_strict_vs_nonstrict_comparative_summary.json",
        "generated/p1233_w1_strict_lane_readiness_gate_summary.json",
    ]

    out = {
        "packet": "P1234",
        "as_of": "2026-05-11",
        "strict_lane_readiness_gate_pass": gate_pass,
        "frozen_artifact_set": artifact_set if gate_pass else [],
        "freeze_note": "Frozen strict-lane artifact baseline for next formalization step." if gate_pass else "Gate not passed; no freeze set declared.",
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1234] gate_pass={gate_pass} wrote {args.out}")


if __name__ == "__main__":
    main()
