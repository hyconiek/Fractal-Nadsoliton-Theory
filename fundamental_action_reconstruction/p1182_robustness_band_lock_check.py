#!/usr/bin/env python3
"""P1182 robustness band lock check for strict-E2E candidate admission."""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    src = GEN / "p1181_multithreshold_delta_curve_summary.json"
    data = json.loads(src.read_text(encoding="utf-8"))

    # Operational lock band from P1181 empirical advantage window.
    band = [0.60, 0.65]
    rows = data.get("rows", [])

    checks = []
    for t in band:
        row = next((r for r in rows if abs(float(r.get("threshold", -999)) - t) < 1e-9), None)
        checks.append({
            "threshold": t,
            "row_found": row is not None,
            "uplift_pass": bool((row or {}).get("uplift", {}).get("integrated_pass")),
            "baseline_pass": bool((row or {}).get("baseline", {}).get("integrated_pass")),
        })

    band_lock_pass = all(c["row_found"] and c["uplift_pass"] and (not c["baseline_pass"]) for c in checks)

    out = {
        "packet": "P1182",
        "as_of": "2026-05-10",
        "source": str(src),
        "band": band,
        "checks": checks,
        "band_lock_pass": band_lock_pass,
        "note": "Band-lock check only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1182_robustness_band_lock_check_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1182] band_lock_pass={band_lock_pass} wrote {out_path}")


if __name__ == "__main__":
    main()
