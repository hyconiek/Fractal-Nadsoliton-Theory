#!/usr/bin/env python3
"""P1237: emit lock report for frozen strict-lane baseline artifacts."""
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
    parser.add_argument("--p1234", type=Path, default=GEN / "p1234_w1_strict_lane_artifact_freeze_note_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1237_w1_frozen_baseline_lock_report_summary.json")
    args = parser.parse_args()

    p1234 = json.loads(args.p1234.read_text(encoding="utf-8"))
    artifacts = p1234.get("frozen_artifact_set", []) if isinstance(p1234.get("frozen_artifact_set", []), list) else []

    lock_rows = []
    for rel in artifacts:
        path = ROOT / rel
        payload = json.loads(path.read_text(encoding="utf-8"))
        lock_rows.append(
            {
                "artifact": rel,
                "packet": payload.get("packet"),
                "theory_closure_status": payload.get("theory_closure_status"),
                "sha256": _sha256(path),
            }
        )

    out = {
        "packet": "P1237",
        "as_of": "2026-05-11",
        "locked_artifacts": lock_rows,
        "lock_count": len(lock_rows),
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Frozen baseline lock report for strict-lane formalization continuity.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1237] lock_count={len(lock_rows)} wrote {args.out}")


if __name__ == "__main__":
    main()
