#!/usr/bin/env python3
"""P1186 repeated post-certification stability rate check."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    repeats = 5
    rows = []
    for i in range(repeats):
        cmd = [sys.executable, str(ROOT / "p1185_post_certification_drift_check.py")]
        p = subprocess.run(cmd, capture_output=True, text=True)
        s = json.loads((GEN / "p1185_post_certification_drift_check_summary.json").read_text(encoding="utf-8"))
        rows.append({
            "iteration": i + 1,
            "returncode": p.returncode,
            "post_certification_stable": bool(s.get("post_certification_stable")),
        })

    stable_count = sum(1 for r in rows if r["returncode"] == 0 and r["post_certification_stable"])
    out = {
        "packet": "P1186",
        "as_of": "2026-05-10",
        "repeats": repeats,
        "stable_count": stable_count,
        "stability_rate": stable_count / repeats if repeats else None,
        "rows": rows,
        "note": "Repeated stability-rate check only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1186_post_certification_stability_rate_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1186] stability_rate={out['stability_rate']} wrote {out_path}")

if __name__ == "__main__":
    main()
