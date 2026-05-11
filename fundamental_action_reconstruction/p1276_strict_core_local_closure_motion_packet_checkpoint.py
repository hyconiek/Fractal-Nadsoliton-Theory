#!/usr/bin/env python3
"""P1276: strict-core local closure motion packet with explicit QW-2191 status."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1275", type=Path, default=GEN / "p1275_strict_kernel_only_scope_and_legacy_historical_only_summary.json")
    parser.add_argument("--p1274", type=Path, default=GEN / "p1274_strict_core_independent_second_pass_symbolic_audit_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1276_strict_core_local_closure_motion_packet_summary.json")
    args = parser.parse_args()

    p1275 = _read(args.p1275)
    p1274 = _read(args.p1274)

    if p1275.get("lane") != "STRICT_CORE_ONLY" or p1274.get("lane") != "STRICT_CORE_ONLY":
        raise SystemExit("P1276 requires STRICT_CORE_ONLY inputs.")

    qw2191_status = "NOT_CLOSED"

    out = {
        "packet": "P1276",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1275": str(args.p1275), "p1274": str(args.p1274)},
        "local_motion": {
            "status": "LOCAL_STRICT_PROOF_READY",
            "scope": "strict-kernel local formal chain",
        },
        "qw2191_closure_status": qw2191_status,
        "global_closure_status": "OPEN",
        "blocking_conditions": [
            "QW-2191 full closure theorem not discharged.",
            "B1/NB1 governance theorem not discharged.",
        ],
        "closure_policy": "NO_GLOBAL_CLOSURE_CLAIM",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1276] wrote {args.out}; qw2191_status={qw2191_status}")


if __name__ == "__main__":
    main()
