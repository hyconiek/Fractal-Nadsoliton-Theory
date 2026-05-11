#!/usr/bin/env python3
"""P1180 compare strict E2E outcomes for baseline vs uplift candidate."""
from __future__ import annotations
import json, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def run_case(candidate: Path, threshold: float) -> dict:
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as tf:
        reg_path = Path(tf.name)
        tf.write(json.dumps({"candidates": [str(candidate)]}))
    cmd = [
        sys.executable, str(ROOT / "p1169_e2e_candidate_integration.py"),
        "--strict-e2e", "--require-out-of-locality-robustness",
        "--robustness-threshold", str(threshold),
        "--candidate", str(candidate),
        "--registry", str(reg_path),
    ]
    p = subprocess.run(cmd, capture_output=True, text=True)
    s = json.loads((GEN / "p1169_e2e_candidate_integration_summary.json").read_text(encoding="utf-8"))
    return {
        "candidate": str(candidate),
        "threshold": threshold,
        "returncode": p.returncode,
        "integrated_pass": bool(s.get("integrated_pass")),
        "robust_fraction": ((s.get("out_of_locality_robustness_summary") or {}).get("robust_fraction")),
    }


def main() -> None:
    threshold = 0.60
    base = (GEN / "p1170_neighbor_1.json").resolve()
    uplift = (GEN / "p1170_neighbor_3.json").resolve()
    rows = [run_case(base, threshold), run_case(uplift, threshold)]
    rows.sort(key=lambda r: (-(r.get("robust_fraction") or -1), r["candidate"]))
    out = {
        "packet": "P1180",
        "as_of": "2026-05-10",
        "threshold": threshold,
        "rows": rows,
        "winner": rows[0] if rows else None,
        "note": "Delta comparison only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1180_e2e_candidate_delta_comparison_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1180] threshold={threshold} winner={out['winner']['candidate'] if out['winner'] else None} wrote {out_path}")

if __name__ == "__main__":
    main()
