#!/usr/bin/env python3
"""P1303: R19 NB1 non-transfer lemma pack v1 checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1302", type=Path, default=GEN / "p1302_qw2191_r18_nb1_nontransfer_obligation_matrix_and_proof_sketch_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1303_qw2191_r19_nb1_nontransfer_lemma_pack_v1_summary.json")
    args = parser.parse_args()

    p1302 = _read(args.p1302)
    if p1302.get("next_priority") != "R19_NB1_NONTRANSFER_LEMMA_PACK_V1":
        raise SystemExit("P1303 requires next_priority=R19_NB1_NONTRANSFER_LEMMA_PACK_V1 from P1302.")

    if p1302.get("proof_sketch_status") != "MATRIX_DRAFTED":
        raise SystemExit("P1303 requires MATRIX_DRAFTED proof sketch status from P1302.")

    lemmas = {
        "LNB1_1": {
            "title": "Role-domain non-isomorphism lemma",
            "maps_to_obligation": "O1_DOMAIN_MISMATCH",
            "status": "DRAFTED",
        },
        "LNB1_2": {
            "title": "No strict-preserving transport functor lemma",
            "maps_to_obligation": "O2_TRANSPORT_FUNCTOR_ABSENCE",
            "status": "DRAFTED",
        },
        "LNB1_3": {
            "title": "QW-2191 persistence under NB assumptions lemma",
            "maps_to_obligation": "O3_QW2191_PERSISTENCE",
            "status": "DRAFTED",
        },
    }

    out = {
        "packet": "P1303",
        "as_of": "2026-05-11",
        "lane": "NB1_NONBRIDGE_TRACK",
        "input": {"p1302": str(args.p1302)},
        "r19_lemma_pack": lemmas,
        "lemma_pack_status": "V1_DRAFT_COMPLETE",
        "next_priority": "R20_NB1_NONTRANSFER_LEMMA_PROOF_ATTEMPT_BATCH_1",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1303] wrote {args.out}; lemmas={len(lemmas)}")


if __name__ == "__main__":
    main()
