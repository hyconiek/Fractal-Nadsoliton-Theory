#!/usr/bin/env python3
"""P1243: semantic consistency check between strict statement and traced assumptions."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _contains_all(text: str, keywords: list[str]) -> bool:
    txt = text.lower()
    return all(k.lower() in txt for k in keywords)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scaffold", type=Path, default=GEN / "p1239_w1_strict_formal_step0_scaffold.json")
    parser.add_argument("--trace", type=Path, default=GEN / "p1242_w1_strict_assumption_trace_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1243_w1_strict_statement_consistency_summary.json")
    args = parser.parse_args()

    scaffold = json.loads(args.scaffold.read_text(encoding="utf-8"))
    trace = json.loads(args.trace.read_text(encoding="utf-8"))

    statement = str(scaffold.get("w1_strict_formal_statement", ""))
    assumptions = scaffold.get("w1_strict_formal_assumptions", []) if isinstance(scaffold.get("w1_strict_formal_assumptions", []), list) else []

    checks = {
        "mentions_strict_lane": _contains_all(statement, ["strict lane"]),
        "mentions_strict_reference": _contains_all(statement, ["strict", "reference"]),
        "mentions_non_strict_exclusion": _contains_all(statement, ["without", "non-strict"]),
        "assumptions_traced": bool(trace.get("all_assumptions_traced", False)),
        "assumption_count_matches": int(trace.get("assumption_count", -1)) == len(assumptions),
    }

    consistent = all(checks.values())

    out = {
        "packet": "P1243",
        "as_of": "2026-05-11",
        "statement_consistency_checks": checks,
        "strict_statement_consistent_with_trace": consistent,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Semantic consistency checkpoint for strict step-0 statement against traced assumptions.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1243] consistent={consistent} wrote {args.out}")


if __name__ == "__main__":
    main()
