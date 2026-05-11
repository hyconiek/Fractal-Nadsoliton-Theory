#!/usr/bin/env python3
"""P1305: R21 NB1 lemma-gap closure provider and invariant packet checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1304", type=Path, default=GEN / "p1304_qw2191_r20_nb1_nontransfer_lemma_proof_attempt_batch_1_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1305_qw2191_r21_nb1_lemma_gap_closure_provider_and_invariant_packet_summary.json")
    args = parser.parse_args()

    p1304 = _read(args.p1304)
    if p1304.get("next_priority") != "R21_NB1_LEMMA_GAP_CLOSURE_PROVIDER_AND_INVARIANT_PACKET":
        raise SystemExit("P1305 requires next_priority=R21_NB1_LEMMA_GAP_CLOSURE_PROVIDER_AND_INVARIANT_PACKET from P1304.")
    if p1304.get("batch_status") != "PARTIAL_PROGRESS":
        raise SystemExit("P1305 requires PARTIAL_PROGRESS status from P1304.")

    packet = {
        "provider_schema_v1": {
            "strict_role_codomain": ["strict_observable", "strict_operator", "strict_selector_token"],
            "legacy_role_codomain": ["legacy_effective_parameter", "historical_fit_token"],
            "non_isomorphism_claim": "No bijection preserves closure predicates across codomains.",
        },
        "strict_invariant_set_v1": [
            "selector_uniqueness_guard",
            "qw2191_obstruction_consistency",
            "no_implicit_legacy_role_transport",
        ],
        "gap_closure_targets": ["LNB1_1", "LNB1_2"],
        "status": "PACKET_READY_FOR_BATCH2",
    }

    out = {
        "packet": "P1305",
        "as_of": "2026-05-11",
        "lane": "NB1_NONBRIDGE_TRACK",
        "input": {"p1304": str(args.p1304)},
        "r21_gap_closure_packet": packet,
        "next_priority": "R22_NB1_NONTRANSFER_LEMMA_PROOF_ATTEMPT_BATCH_2",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1305] wrote {args.out}; status={packet['status']}")


if __name__ == "__main__":
    main()
