#!/usr/bin/env python3
"""P1272: strict-core epsilon formal-proof text and symbolic-verifier checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1271", type=Path, default=GEN / "p1271_strict_core_epsilon_derivation_theorem_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1272_strict_core_epsilon_formal_proof_and_symbolic_verifier_summary.json")
    args = parser.parse_args()

    p1271 = _read(args.p1271)
    if p1271.get("lane") != "STRICT_CORE_ONLY":
        raise SystemExit("P1272 requires STRICT_CORE_ONLY lane from P1271.")

    proof_text_status = "DRAFT_V1"
    symbolic_verifier_status = "PARTIAL_RUN"

    out = {
        "packet": "P1272",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1271": str(args.p1271)},
        "formal_proof_text": {
            "theorem_id": "THM_EPS_QW2191",
            "version": "v1",
            "status": proof_text_status,
            "coverage": ["D1", "D2", "D3"],
        },
        "symbolic_verifier": {
            "id": "SV_EPS_QW2191",
            "status": symbolic_verifier_status,
            "verified_steps": ["D1", "D3"],
            "pending_steps": ["D2"],
        },
        "remaining_obligations": [
            "Complete symbolic verification of D2 decomposition step.",
            "Promote proof text from DRAFT_V1 to REVIEWED_FORMAL_PROOF.",
        ],
        "strict_kernel_closure_ready": False,
        "closure_policy": "STRICT_KERNEL_CLOSURE_FORBIDDEN_UNTIL_SYMBOLIC_VERIFIER_COMPLETE_AND_B1_OR_NB1",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1272] wrote {args.out}; proof={proof_text_status} verifier={symbolic_verifier_status}")


if __name__ == "__main__":
    main()
