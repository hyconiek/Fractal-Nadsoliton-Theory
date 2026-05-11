#!/usr/bin/env python3
"""P1240: check whether strict formal step-0 scaffold placeholders are filled."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _contains_todo(value: object) -> bool:
    if isinstance(value, str):
        return "TODO_" in value
    if isinstance(value, list):
        return any(_contains_todo(v) for v in value)
    if isinstance(value, dict):
        return any(_contains_todo(v) for v in value.values())
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scaffold", type=Path, default=GEN / "p1239_w1_strict_formal_step0_scaffold.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1240_w1_strict_statement_fill_summary.json")
    args = parser.parse_args()

    scaffold = json.loads(args.scaffold.read_text(encoding="utf-8"))

    statement = scaffold.get("w1_strict_formal_statement", "")
    assumptions = scaffold.get("w1_strict_formal_assumptions", [])

    statement_filled = isinstance(statement, str) and statement != "" and "TODO_" not in statement
    assumptions_filled = isinstance(assumptions, list) and len(assumptions) > 0 and not _contains_todo(assumptions)
    placeholders_present = _contains_todo({"statement": statement, "assumptions": assumptions})

    out = {
        "packet": "P1240",
        "as_of": "2026-05-11",
        "statement_filled": statement_filled,
        "assumptions_filled": assumptions_filled,
        "placeholders_present": placeholders_present,
        "strict_formal_step0_ready": statement_filled and assumptions_filled and not placeholders_present,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Fill-check checkpoint for strict formal step-0 scaffold.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1240] ready={out['strict_formal_step0_ready']} wrote {args.out}")


if __name__ == "__main__":
    main()
