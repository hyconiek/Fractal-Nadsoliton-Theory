#!/usr/bin/env python3
"""P1285: R4 strict selector source identification checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1284", type=Path, default=GEN / "p1284_qw2191_r3_independent_audit_and_replication_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1285_qw2191_r4_strict_selector_source_identification_summary.json")
    args = parser.parse_args()

    p1284 = _read(args.p1284)
    if p1284.get("next_priority") != "R4_STRICT_SELECTOR_SOURCE_IDENTIFICATION":
        raise SystemExit("P1285 requires next_priority=R4_STRICT_SELECTOR_SOURCE_IDENTIFICATION from P1284.")
    if p1284.get("independent_audit", {}).get("result") != "PASS":
        raise SystemExit("P1285 requires PASS independent audit from P1284.")

    source_classes = [
        {
            "id": "SSEL_SRC_A",
            "description": "asymmetric strict-boundary phase source",
            "identifiability": "DISTINGUISHABLE",
        },
        {
            "id": "SSEL_SRC_B",
            "description": "strict-lane transport-induced selector seed",
            "identifiability": "DISTINGUISHABLE",
        },
    ]

    out = {
        "packet": "P1285",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1284": str(args.p1284)},
        "selector_source_identification": {
            "candidate_classes": source_classes,
            "minimal_source_family": ["SSEL_SRC_A", "SSEL_SRC_B"],
            "nonuniqueness_residual": "OPEN",
            "status": "PARTIAL_DISCHARGE",
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R5_SELECTOR_NONUNIQUENESS_EXCLUSION_ATTEMPT",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1285] wrote {args.out}; classes={len(source_classes)} status=PARTIAL_DISCHARGE")


if __name__ == "__main__":
    main()
