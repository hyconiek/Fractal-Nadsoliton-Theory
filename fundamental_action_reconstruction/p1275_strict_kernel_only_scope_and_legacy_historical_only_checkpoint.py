#!/usr/bin/env python3
"""P1275: strict-kernel-only scope declaration with legacy set to historical context only."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1274", type=Path, default=GEN / "p1274_strict_core_independent_second_pass_symbolic_audit_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1275_strict_kernel_only_scope_and_legacy_historical_only_summary.json")
    args = parser.parse_args()

    p1274 = _read(args.p1274)
    if p1274.get("lane") != "STRICT_CORE_ONLY":
        raise SystemExit("P1275 requires STRICT_CORE_ONLY lane from P1274.")

    out = {
        "packet": "P1275",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1274": str(args.p1274)},
        "scope_policy": {
            "active_kernel": "K_strict_gate",
            "legacy_kernel_role": "HISTORICAL_NADSOLITON_DEFINITION_ATTEMPT_ONLY",
            "legacy_kernel_operational_use": "DISALLOWED_IN_CURRENT_STRICT_CORE_PROOF_CHAIN",
        },
        "methodological_constraints": [
            "No legacy-role transfer into strict-core proof steps.",
            "All closure-relevant lemmas must be strict-kernel internal.",
            "Global closure still requires explicit B1/NB1 governance decision.",
        ],
        "strict_kernel_closure_ready": False,
        "closure_policy": "STRICT_ONLY_FORMAL_CHAIN_CONTINUES__GLOBAL_CLOSURE_GATED",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1275] wrote {args.out}; active_kernel=K_strict_gate")


if __name__ == "__main__":
    main()
