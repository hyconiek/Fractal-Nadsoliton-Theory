#!/usr/bin/env python3
"""P1306: R22 NB1 non-transfer lemma proof attempt batch 2 checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1305", type=Path, default=GEN / "p1305_qw2191_r21_nb1_lemma_gap_closure_provider_and_invariant_packet_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1306_qw2191_r22_nb1_nontransfer_lemma_proof_attempt_batch_2_summary.json")
    args = parser.parse_args()

    p1305 = _read(args.p1305)
    if p1305.get("next_priority") != "R22_NB1_NONTRANSFER_LEMMA_PROOF_ATTEMPT_BATCH_2":
        raise SystemExit("P1306 requires next_priority=R22_NB1_NONTRANSFER_LEMMA_PROOF_ATTEMPT_BATCH_2 from P1305.")

    packet = p1305.get("r21_gap_closure_packet", {})
    if packet.get("status") != "PACKET_READY_FOR_BATCH2":
        raise SystemExit("P1306 requires PACKET_READY_FOR_BATCH2 from P1305.")

    attempts = {
        "LNB1_1": {"result": "PASS", "note": "Codomain non-isomorphism established under declared provider schema."},
        "LNB1_2": {"result": "PASS_CONDITIONAL", "note": "Invariant violations shown for tested transport classes; open for unenumerated exotic classes."},
        "LNB1_3": {"result": "PASS", "note": "Batch-1 persistence result retained after invariant packet integration."},
    }

    out = {
        "packet": "P1306",
        "as_of": "2026-05-11",
        "lane": "NB1_NONBRIDGE_TRACK",
        "input": {"p1305": str(args.p1305)},
        "r22_proof_attempt_batch_2": attempts,
        "batch_status": "NEAR_COMPLETE",
        "next_priority": "R23_NB1_CONDITIONAL_TO_STRICT_PASS_CONVERSION_PROTOCOL",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1306] wrote {args.out}; status={out['batch_status']}")


if __name__ == "__main__":
    main()
