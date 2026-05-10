#!/usr/bin/env python3
"""P1193: First strict-rigor discharge attempt for W4 defect-zero witness."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    p1192 = json.loads((GEN / "p1192_selector_premise_witness_map_summary.json").read_text(encoding="utf-8"))

    # Conservative protocol constants (no false pass):
    tol_strict = 1e-10

    # Current repo state has no exported dedicated W4 symbolic-zero certificate;
    # therefore this attempt remains bounded and negative by default.
    observed_max_abs_defect = None
    has_symbolic_certificate = False

    w4_pass = bool(has_symbolic_certificate and observed_max_abs_defect is not None and observed_max_abs_defect <= tol_strict)

    out = {
        "packet": "P1193",
        "as_of": "2026-05-10",
        "input_open_count": int(p1192.get("open_count", -1)),
        "target_witness": "W4_defect_polynomial_zeroes",
        "tolerance_strict": tol_strict,
        "has_symbolic_certificate": has_symbolic_certificate,
        "observed_max_abs_defect": observed_max_abs_defect,
        "w4_discharge_pass": w4_pass,
        "status": "OPEN" if not w4_pass else "DISCHARGED",
        "note": "No exported symbolic-zero certificate yet; strict no-false-pass keeps W4 open.",
    }

    out_path = GEN / "p1193_w4_defect_zero_witness_attempt_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1193] w4_discharge_pass={w4_pass} status={out['status']} wrote {out_path}")


if __name__ == "__main__":
    main()
