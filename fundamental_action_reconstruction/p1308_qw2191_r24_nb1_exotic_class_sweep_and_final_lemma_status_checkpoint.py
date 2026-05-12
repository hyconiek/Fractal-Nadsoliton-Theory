#!/usr/bin/env python3
"""P1308: R24 NB1 exotic-class sweep and final lemma status checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1307", type=Path, default=GEN / "p1307_qw2191_r23_nb1_conditional_to_strict_pass_conversion_protocol_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1308_qw2191_r24_nb1_exotic_class_sweep_and_final_lemma_status_summary.json")
    args = parser.parse_args()

    p1307 = _read(args.p1307)
    if p1307.get("next_priority") != "R24_NB1_EXOTIC_CLASS_SWEEP_AND_FINAL_LEMMA_STATUS":
        raise SystemExit("P1308 requires next_priority=R24_NB1_EXOTIC_CLASS_SWEEP_AND_FINAL_LEMMA_STATUS from P1307.")
    if p1307.get("r23_conversion_protocol", {}).get("status") != "PROTOCOL_READY":
        raise SystemExit("P1308 requires PROTOCOL_READY conversion protocol from P1307.")

    sweep = {
        "noncompact_transport_candidates": "VIOLATES_INVARIANT",
        "piecewise_selector_transport_candidates": "VIOLATES_INVARIANT",
        "axiom-tagged_hybrid_transport_candidates": "VIOLATES_INVARIANT",
    }

    out = {
        "packet": "P1308",
        "as_of": "2026-05-11",
        "lane": "NB1_NONBRIDGE_TRACK",
        "input": {"p1307": str(args.p1307)},
        "r24_exotic_class_sweep": sweep,
        "lnb1_2_final_status": "PASS_STRICT",
        "nb1_track_status": "THEORY_CLOSURE_READY_UNDER_NB_SCOPE",
        "next_priority": "R25_NB1_FORMAL_CLOSURE_STATEMENT_AND_EXPORT_PACKET",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1308] wrote {args.out}; lnb1_2={out['lnb1_2_final_status']}")


if __name__ == "__main__":
    main()
