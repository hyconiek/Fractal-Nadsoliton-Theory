#!/usr/bin/env python3
"""P1280: R2 bound-transport checkpoint after P1279 completeness status."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1279", type=Path, default=GEN / "p1279_qw2191_r1_completeness_theorem_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1280_qw2191_r2_bound_transport_checkpoint_summary.json")
    args = parser.parse_args()

    p1279 = _read(args.p1279)
    if p1279.get("theorem", {}).get("status") != "PARTIAL_DISCHARGE":
        raise SystemExit("P1280 requires P1279 theorem status PARTIAL_DISCHARGE.")
    if p1279.get("next_priority") != "R2_BOUND_TRANSPORT":
        raise SystemExit("P1280 requires P1279 next_priority=R2_BOUND_TRANSPORT.")

    obligations = [
        {
            "id": "R2.O1",
            "name": "Transport all declared QW-2191 local bounds to one common gauge.",
            "status": "OPEN",
        },
        {
            "id": "R2.O2",
            "name": "Provide mismatch-control lemma for transported bounds.",
            "status": "OPEN",
        },
        {
            "id": "R2.O3",
            "name": "Export machine-checkable transport certificate.",
            "status": "OPEN",
        },
    ]

    out = {
        "packet": "P1280",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1279": str(args.p1279)},
        "r2_program": {
            "name": "QW2191_BOUND_TRANSPORT",
            "entry_gate": "OPEN",
            "obligations": obligations,
            "all_obligations_discharged": False,
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
            "bridge_decision_required": True,
        },
        "next_priority": "R2_O1_COMMON_GAUGE_TRANSPORT",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1280] wrote {args.out}; obligations_open={len(obligations)}")


if __name__ == "__main__":
    main()
