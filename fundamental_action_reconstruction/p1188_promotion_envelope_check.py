#!/usr/bin/env python3
"""P1188 promotion envelope check from stress-profile matrix."""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    src = GEN / "p1187_stability_stress_profile_matrix_summary.json"
    m = json.loads(src.read_text(encoding="utf-8"))

    uplift_key = "uplift_neighbor_3"
    baseline_key = "baseline_neighbor_1"
    required_band = [0.60, 0.65]
    required_rate = 0.95

    matrix = m.get("matrix", {})
    uplift = matrix.get(uplift_key, {})
    baseline = matrix.get(baseline_key, {})

    checks = []
    for t in required_band:
        tk = str(t)
        ur = float((uplift.get(tk) or {}).get("stability_rate", -1.0))
        br = float((baseline.get(tk) or {}).get("stability_rate", -1.0))
        checks.append({
            "threshold": t,
            "uplift_rate": ur,
            "baseline_rate": br,
            "uplift_meets": ur >= required_rate,
            "baseline_fails": br < required_rate,
        })

    promote = all(c["uplift_meets"] and c["baseline_fails"] for c in checks)

    out = {
        "packet": "P1188",
        "as_of": "2026-05-10",
        "source": str(src),
        "uplift_key": uplift_key,
        "baseline_key": baseline_key,
        "required_band": required_band,
        "required_rate": required_rate,
        "checks": checks,
        "promotion_envelope_pass": promote,
        "note": "Envelope check only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1188_promotion_envelope_check_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1188] promotion_envelope_pass={promote} wrote {out_path}")

if __name__ == "__main__":
    main()
