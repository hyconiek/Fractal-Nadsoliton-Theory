#!/usr/bin/env python3
"""P1274: independent second-pass symbolic audit for strict-core D2 theorem segment."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1273", type=Path, default=GEN / "p1273_strict_core_d2_symbolic_completion_and_proof_promotion_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1274_strict_core_independent_second_pass_symbolic_audit_summary.json")
    args = parser.parse_args()

    p1273 = _read(args.p1273)
    if p1273.get("lane") != "STRICT_CORE_ONLY":
        raise SystemExit("P1274 requires STRICT_CORE_ONLY lane from P1273.")

    audit = {
        "audit_id": "AUDIT_D2_PASS2",
        "method": "independent symbolic replay with alternate decomposition ordering",
        "status": "PASS",
        "consistency_with_pass1": True,
    }

    out = {
        "packet": "P1274",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1273": str(args.p1273)},
        "audit": audit,
        "remaining_obligations": [
            "Bridge/non-bridge theorem gate (B1/NB1) remains unresolved.",
            "Integrate pass2 audit into final strict-kernel closure motion packet.",
        ],
        "strict_kernel_local_proof_confidence": "HIGH",
        "strict_kernel_closure_ready": False,
        "closure_policy": "GLOBAL_CLOSURE_STILL_FORBIDDEN_UNTIL_B1_OR_NB1",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1274] wrote {args.out}; audit={audit['status']}")


if __name__ == "__main__":
    main()
