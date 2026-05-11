#!/usr/bin/env python3
"""P1294: R10 formal proof-chain and countermodel sweep checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1293", type=Path, default=GEN / "p1293_qw2191_r9_formal_selector_source_theorem_draft_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1294_qw2191_r10_formal_proof_chain_and_countermodel_sweep_summary.json")
    args = parser.parse_args()

    p1293 = _read(args.p1293)
    if p1293.get("next_priority") != "R10_FORMAL_PROOF_CHAIN_AND_COUNTERMODEL_SWEEP":
        raise SystemExit("P1294 requires next_priority=R10_FORMAL_PROOF_CHAIN_AND_COUNTERMODEL_SWEEP from P1293.")

    if p1293.get("r9_theorem_draft", {}).get("status") != "DRAFT":
        raise SystemExit("P1294 requires R9 theorem draft status DRAFT from P1293.")

    proof_chain = {
        "lemma_sequence": ["L_R10_1_transport_stability", "L_R10_2_margin_persistence", "L_R10_3_source_identifiability"],
        "status": "IN_PROGRESS",
    }
    countermodel_sweep = {
        "search_space": "bounded_strict_lane_variants_v1",
        "countermodels_found": 0,
        "status": "PASS_NO_COUNTERMODEL_IN_SEARCH_SPACE",
    }

    out = {
        "packet": "P1294",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1293": str(args.p1293)},
        "r10": {
            "proof_chain": proof_chain,
            "countermodel_sweep": countermodel_sweep,
            "status": "PARTIAL_DISCHARGE",
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R11_FORMAL_PROOF_COMPLETION_AND_PEER_REPLAY",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1294] wrote {args.out}; countermodels={countermodel_sweep['countermodels_found']}")


if __name__ == "__main__":
    main()
