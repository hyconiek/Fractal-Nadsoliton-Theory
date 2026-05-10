#!/usr/bin/env python3
"""P1190 strict nonclosure guard for promotion artifacts."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

BANNED_TOKENS = (
    "toe closure achieved",
    "theory closed",
    "qw-2191 discharged",
    "strict-core selector closure achieved",
)


def main() -> None:
    p1189 = json.loads((GEN / "p1189_final_promotion_contract_summary.json").read_text(encoding="utf-8"))

    contract_pass = bool(p1189.get("final_promotion_contract_pass", False))
    note = str(p1189.get("note", "")).strip()
    note_lower = note.lower()

    has_nonclosure_marker = ("no" in note_lower and "closure" in note_lower) or ("qw-2191" in note_lower)
    banned_hits = [token for token in BANNED_TOKENS if token in note_lower]

    guard_pass = contract_pass and has_nonclosure_marker and not banned_hits

    out = {
        "packet": "P1190",
        "as_of": "2026-05-10",
        "input_contract_pass": contract_pass,
        "has_explicit_nonclosure_marker": has_nonclosure_marker,
        "banned_claim_hits": banned_hits,
        "strict_promotion_nonclosure_guard_pass": guard_pass,
        "note": "Operational promotion may pass, but closure claims remain forbidden.",
    }

    out_path = GEN / "p1190_strict_promotion_nonclosure_guard_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1190] strict_promotion_nonclosure_guard_pass={guard_pass} wrote {out_path}")


if __name__ == "__main__":
    main()
