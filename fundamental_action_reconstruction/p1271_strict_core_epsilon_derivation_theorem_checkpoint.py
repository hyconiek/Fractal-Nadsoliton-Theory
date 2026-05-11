#!/usr/bin/env python3
"""P1271: theorem-grade derivation checkpoint for epsilon_qw2191 in strict-core."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1270", type=Path, default=GEN / "p1270_strict_core_qw2191_full_obstruction_bound_proof_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1271_strict_core_epsilon_derivation_theorem_summary.json")
    args = parser.parse_args()

    p1270 = _read(args.p1270)
    if p1270.get("lane") != "STRICT_CORE_ONLY":
        raise SystemExit("P1271 requires STRICT_CORE_ONLY lane from P1270.")

    theorem = {
        "id": "THM_EPS_QW2191",
        "statement": "epsilon_qw2191 bound is derivable from strict-core primitives with no external selector premise.",
        "status": "PARTIAL_DISCHARGE",
    }

    derivation_chain = [
        "D1: strict-core locality envelope controls residual growth.",
        "D2: obstruction interface decomposition yields bounded residual contribution.",
        "D3: assembled bound implies epsilon_qw2191 <= 0.047 on declared domain.",
    ]

    unresolved = [
        "Need full formal proof text replacing derivation sketch D1-D3.",
        "Need independent symbolic verification artifact for D2 decomposition.",
    ]

    out = {
        "packet": "P1271",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1270": str(args.p1270)},
        "theorem": theorem,
        "derivation_chain": derivation_chain,
        "unresolved": unresolved,
        "strict_kernel_closure_ready": False,
        "closure_policy": "STRICT_KERNEL_CLOSURE_FORBIDDEN_UNTIL_THM_EPS_QW2191_DISCHARGED_AND_B1_OR_NB1",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1271] wrote {args.out}; status={theorem['status']}")


if __name__ == "__main__":
    main()
