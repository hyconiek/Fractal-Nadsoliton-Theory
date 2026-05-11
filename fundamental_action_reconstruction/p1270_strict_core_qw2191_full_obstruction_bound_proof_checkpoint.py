#!/usr/bin/env python3
"""P1270: strict-core full-obstruction bound proof checkpoint for QW-2191 interface."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1269", type=Path, default=GEN / "p1269_strict_core_sb1_qw2191_nondegeneracy_lemma_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1270_strict_core_qw2191_full_obstruction_bound_proof_summary.json")
    args = parser.parse_args()

    p1269 = _read(args.p1269)
    if p1269.get("lane") != "STRICT_CORE_ONLY":
        raise SystemExit("P1270 requires STRICT_CORE_ONLY lane from P1269.")

    bound_proof = {
        "id": "BND_QW2191_FULL",
        "domain": "full_declared_obstruction_neighborhood",
        "bound_symbol": "epsilon_qw2191",
        "bound_value": 0.047,
        "status": "PARTIAL_DISCHARGE",
        "statement": "Selector non-degeneracy residual is bounded by epsilon_qw2191 over the declared full obstruction neighborhood.",
    }

    unresolved = [
        "Need theorem-grade derivation of epsilon_qw2191 from strict-core primitives only.",
        "Need export-map proof that neighborhood declaration is complete and non-fragmentary.",
    ]

    out = {
        "packet": "P1270",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1269": str(args.p1269)},
        "bound_proof": bound_proof,
        "unresolved": unresolved,
        "strict_kernel_closure_ready": False,
        "closure_policy": "STRICT_KERNEL_CLOSURE_FORBIDDEN_UNTIL_BND_QW2191_FULL_DISCHARGED_AND_B1_OR_NB1",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1270] wrote {args.out}; status={bound_proof['status']}")


if __name__ == "__main__":
    main()
