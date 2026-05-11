#!/usr/bin/env python3
"""P1179 targeted out-of-locality uplift scan over neighbor candidates."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def run_probe(candidate_path: Path) -> dict:
    cmd = [sys.executable, str(ROOT / "p1171_out_of_locality_robustness_probe.py"), str(candidate_path)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    s = json.loads((GEN / "p1171_out_of_locality_robustness_probe_summary.json").read_text(encoding="utf-8"))
    return {
        "candidate": str(candidate_path),
        "returncode": p.returncode,
        "robust_fraction": s.get("robust_fraction"),
        "robust_cases": s.get("robust_cases"),
        "cases": s.get("cases"),
    }


def main() -> None:
    candidates = [GEN / f"p1170_neighbor_{i}.json" for i in range(1, 6)]
    rows = [run_probe(c.resolve()) for c in candidates if c.exists()]
    rows.sort(key=lambda r: (-(r.get("robust_fraction") or -1), r["candidate"]))
    best = rows[0] if rows else None
    out = {
        "packet": "P1179",
        "as_of": "2026-05-10",
        "candidates_scanned": len(rows),
        "rows": rows,
        "best_candidate": best,
        "note": "Targeted uplift scan only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1179_targeted_out_of_locality_uplift_scan_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1179] candidates={len(rows)} best={best['candidate'] if best else None} rf={best['robust_fraction'] if best else None} wrote {out_path}")


if __name__ == "__main__":
    main()
