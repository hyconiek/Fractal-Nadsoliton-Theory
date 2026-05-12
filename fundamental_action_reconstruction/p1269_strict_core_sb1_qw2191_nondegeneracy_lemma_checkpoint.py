#!/usr/bin/env python3
"""P1269: strict-core non-degeneracy lemma checkpoint for SB1 under QW-2191."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1268", type=Path, default=GEN / "p1268_strict_core_sb1_qw2191_compatibility_theorem_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1269_strict_core_sb1_qw2191_nondegeneracy_lemma_summary.json")
    args = parser.parse_args()

    p1268 = _read(args.p1268)
    status = p1268.get("theorem", {}).get("status")
    if status not in {"PARTIAL_COMPATIBLE", "DISCHARGED"}:
        raise SystemExit("P1269 requires at least PARTIAL_COMPATIBLE status from P1268.")

    lemma = {
        "id": "L_SB1_QW2191_ND",
        "statement": (
            "Under strict-core constraints, selector choice remains non-degenerate across the declared "
            "QW-2191 obstruction interface neighborhood."
        ),
        "status": "PARTIAL_DISCHARGE",
    }

    obligations = [
        "Provide explicit bound proof for full obstruction neighborhood (not only sampled slices).",
        "Upgrade compatibility chain to theorem-grade derivation with no external selector premise.",
    ]

    out = {
        "packet": "P1269",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1268": str(args.p1268)},
        "lemma": lemma,
        "evidence_anchor": "generated/p1268_strict_core_sb1_qw2191_compatibility_theorem_summary.json",
        "remaining_obligations": obligations,
        "strict_kernel_closure_ready": False,
        "closure_policy": "STRICT_KERNEL_CLOSURE_FORBIDDEN_UNTIL_L_SB1_QW2191_ND_DISCHARGED_AND_B1_OR_NB1",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1269] wrote {args.out}; status={lemma['status']}")


if __name__ == "__main__":
    main()
