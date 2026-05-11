#!/usr/bin/env python3
"""P1304: R20 NB1 non-transfer lemma proof attempt batch 1 checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1303", type=Path, default=GEN / "p1303_qw2191_r19_nb1_nontransfer_lemma_pack_v1_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1304_qw2191_r20_nb1_nontransfer_lemma_proof_attempt_batch_1_summary.json")
    args = parser.parse_args()

    p1303 = _read(args.p1303)
    if p1303.get("next_priority") != "R20_NB1_NONTRANSFER_LEMMA_PROOF_ATTEMPT_BATCH_1":
        raise SystemExit("P1304 requires next_priority=R20_NB1_NONTRANSFER_LEMMA_PROOF_ATTEMPT_BATCH_1 from P1303.")
    if p1303.get("lemma_pack_status") != "V1_DRAFT_COMPLETE":
        raise SystemExit("P1304 requires V1_DRAFT_COMPLETE lemma pack status from P1303.")

    attempts = {
        "LNB1_1": {"result": "OPEN", "note": "Need explicit codomain typing schema for strict-role layer."},
        "LNB1_2": {"result": "OPEN", "note": "Need closure-invariant list with formal violation witnesses."},
        "LNB1_3": {"result": "PASS_CONDITIONAL", "note": "QW-2191 replay obstruction preserved under null-provider assumption."},
    }

    out = {
        "packet": "P1304",
        "as_of": "2026-05-11",
        "lane": "NB1_NONBRIDGE_TRACK",
        "input": {"p1303": str(args.p1303)},
        "r20_proof_attempt_batch_1": attempts,
        "batch_status": "PARTIAL_PROGRESS",
        "next_priority": "R21_NB1_LEMMA_GAP_CLOSURE_PROVIDER_AND_INVARIANT_PACKET",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1304] wrote {args.out}; status={out['batch_status']}")


if __name__ == "__main__":
    main()
