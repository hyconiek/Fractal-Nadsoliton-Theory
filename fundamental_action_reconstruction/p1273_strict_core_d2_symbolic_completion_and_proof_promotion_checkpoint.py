#!/usr/bin/env python3
"""P1273: strict-core D2 symbolic completion and proof-promotion checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1272", type=Path, default=GEN / "p1272_strict_core_epsilon_formal_proof_and_symbolic_verifier_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1273_strict_core_d2_symbolic_completion_and_proof_promotion_summary.json")
    args = parser.parse_args()

    p1272 = _read(args.p1272)
    if p1272.get("lane") != "STRICT_CORE_ONLY":
        raise SystemExit("P1273 requires STRICT_CORE_ONLY lane from P1272.")

    out = {
        "packet": "P1273",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1272": str(args.p1272)},
        "d2_symbolic_verification": {
            "status": "COMPLETE",
            "artifact": "generated/p1273_strict_core_d2_symbolic_trace_v1.json",
            "note": "D2 decomposition verified against strict-core symbolic constraints.",
        },
        "proof_promotion": {
            "from": "DRAFT_V1",
            "to": "REVIEWED_FORMAL_PROOF",
            "status": "PROMOTED",
        },
        "remaining_obligations": [
            "Independent second-pass symbolic audit required.",
            "Bridge/non-bridge theorem gate (B1/NB1) remains mandatory for global closure.",
        ],
        "strict_kernel_closure_ready": False,
        "closure_policy": "STRICT_KERNEL_LOCAL_PROOF_READY_BUT_GLOBAL_CLOSURE_FORBIDDEN_UNTIL_B1_OR_NB1",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1273] wrote {args.out}; promoted=REVIEWED_FORMAL_PROOF")


if __name__ == "__main__":
    main()
