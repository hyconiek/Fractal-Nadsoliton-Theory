#!/usr/bin/env python3
"""P1235: validate frozen strict-lane artifact set consistency."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1234", type=Path, default=GEN / "p1234_w1_strict_lane_artifact_freeze_note_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1235_w1_frozen_set_consistency_summary.json")
    args = parser.parse_args()

    p1234 = _load(args.p1234)
    artifacts = p1234.get("frozen_artifact_set", []) if isinstance(p1234.get("frozen_artifact_set", []), list) else []

    checks = []
    all_ok = True

    for rel in artifacts:
        path = ROOT / rel.replace("generated/", "generated/")
        exists = path.exists()
        packet_ok = False
        closure_field_ok = False

        if exists:
            try:
                payload = _load(path)
                packet_ok = isinstance(payload.get("packet"), str) and len(payload.get("packet")) > 0
                closure_field_ok = "theory_closure_status" in payload
            except Exception:
                packet_ok = False
                closure_field_ok = False

        ok = exists and packet_ok and closure_field_ok
        all_ok = all_ok and ok
        checks.append({
            "artifact": rel,
            "exists": exists,
            "packet_ok": packet_ok,
            "theory_closure_status_present": closure_field_ok,
            "ok": ok,
        })

    out = {
        "packet": "P1235",
        "as_of": "2026-05-11",
        "frozen_set_size": len(artifacts),
        "all_artifacts_consistent": all_ok,
        "checks": checks,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Consistency check over frozen strict-lane artifact set.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1235] all_artifacts_consistent={all_ok} wrote {args.out}")


if __name__ == "__main__":
    main()
