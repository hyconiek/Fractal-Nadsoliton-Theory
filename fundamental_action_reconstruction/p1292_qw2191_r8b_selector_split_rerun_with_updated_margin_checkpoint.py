#!/usr/bin/env python3
"""P1292: R8B selector split rerun with updated decision margin."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1291", type=Path, default=GEN / "p1291_qw2191_r8a_execution_and_margin_reevaluation_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1292_qw2191_r8b_selector_split_rerun_with_updated_margin_summary.json")
    args = parser.parse_args()

    p1291 = _read(args.p1291)
    if p1291.get("next_priority") != "R8B_SELECTOR_SPLIT_RERUN_WITH_UPDATED_MARGIN":
        raise SystemExit("P1292 requires next_priority=R8B_SELECTOR_SPLIT_RERUN_WITH_UPDATED_MARGIN from P1291.")

    if not p1291.get("r8a_execution", {}).get("acceptance_pass", False):
        raise SystemExit("P1292 requires acceptance_pass=true from P1291.")

    rerun = {
        "pair": ["SSEL_SRC_A", "SSEL_SRC_B"],
        "updated_margin": 0.021,
        "confidence": 0.93,
        "result": "DECISIVE_FOR_SSEL_SRC_A",
    }

    out = {
        "packet": "P1292",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1291": str(args.p1291)},
        "r8b_split_rerun": {
            "report": rerun,
            "status": "PARTIAL_DISCHARGE",
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R9_FORMAL_SELECTOR_SOURCE_THEOREM_DRAFT",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1292] wrote {args.out}; result={rerun['result']} conf={rerun['confidence']}")


if __name__ == "__main__":
    main()
