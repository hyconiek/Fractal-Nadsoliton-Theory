#!/usr/bin/env python3
"""P1238: strict formalization entry gate using frozen baseline lock report."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1237", type=Path, default=GEN / "p1237_w1_frozen_baseline_lock_report_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1238_w1_strict_formalization_entry_gate_summary.json")
    args = parser.parse_args()

    lock = json.loads(args.p1237.read_text(encoding="utf-8"))
    rows = lock.get("locked_artifacts", []) if isinstance(lock.get("locked_artifacts", []), list) else []

    checks = []
    all_match = True
    for row in rows:
        rel = row.get("artifact")
        expected = row.get("sha256")
        path = ROOT / rel
        exists = path.exists()
        actual = _sha256(path) if exists else None
        match = exists and actual == expected
        all_match = all_match and match
        checks.append({
            "artifact": rel,
            "exists": exists,
            "expected_sha256": expected,
            "actual_sha256": actual,
            "match": match,
        })

    out = {
        "packet": "P1238",
        "as_of": "2026-05-11",
        "lock_count": len(rows),
        "all_hashes_match": all_match,
        "strict_formalization_entry_gate_pass": all_match,
        "checks": checks,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Entry gate pass means frozen baseline unchanged since P1237 lock report.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1238] gate_pass={all_match} wrote {args.out}")


if __name__ == "__main__":
    main()
