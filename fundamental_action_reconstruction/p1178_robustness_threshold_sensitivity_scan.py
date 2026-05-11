#!/usr/bin/env python3
"""P1178 robustness-threshold sensitivity scan for P1169 strict E2E."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def run(threshold: float) -> dict:
    cmd = [sys.executable, str(ROOT / "p1169_e2e_candidate_integration.py"), "--strict-e2e", "--require-out-of-locality-robustness", "--robustness-threshold", str(threshold)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    s = json.loads((GEN / "p1169_e2e_candidate_integration_summary.json").read_text(encoding="utf-8"))
    return {
        "threshold": threshold,
        "returncode": p.returncode,
        "integrated_pass": bool(s.get("integrated_pass")),
        "robust_fraction": ((s.get("out_of_locality_robustness_summary") or {}).get("robust_fraction")),
    }


def main() -> None:
    thresholds = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99]
    rows = [run(t) for t in thresholds]
    passing = [r for r in rows if r["integrated_pass"]]
    min_pass = min((r["threshold"] for r in passing), default=None)
    out = {
        "packet": "P1178",
        "as_of": "2026-05-10",
        "thresholds": thresholds,
        "rows": rows,
        "min_threshold_with_pass": min_pass,
        "note": "Sensitivity scan only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1178_robustness_threshold_sensitivity_scan_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1178] rows={len(rows)} min_threshold_with_pass={min_pass} wrote {out_path}")


if __name__ == "__main__":
    main()
