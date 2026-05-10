#!/usr/bin/env python3
"""P1181 multi-threshold delta curve for baseline vs uplift candidates."""
from __future__ import annotations
import json, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def run_case(candidate: Path, threshold: float) -> dict:
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as tf:
        reg_path = Path(tf.name)
        tf.write(json.dumps({"candidates": [str(candidate)]}))
    cmd = [sys.executable, str(ROOT / "p1169_e2e_candidate_integration.py"), "--strict-e2e", "--require-out-of-locality-robustness", "--robustness-threshold", str(threshold), "--candidate", str(candidate), "--registry", str(reg_path)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    s = json.loads((GEN / "p1169_e2e_candidate_integration_summary.json").read_text(encoding="utf-8"))
    return {
        "threshold": threshold,
        "candidate": str(candidate),
        "returncode": p.returncode,
        "integrated_pass": bool(s.get("integrated_pass")),
        "robust_fraction": ((s.get("out_of_locality_robustness_summary") or {}).get("robust_fraction")),
    }


def main() -> None:
    thresholds = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75]
    baseline = (GEN / "p1170_neighbor_1.json").resolve()
    uplift = (GEN / "p1170_neighbor_3.json").resolve()

    rows = []
    for t in thresholds:
        b = run_case(baseline, t)
        u = run_case(uplift, t)
        rows.append({"threshold": t, "baseline": b, "uplift": u, "uplift_advantage": bool(u["integrated_pass"]) and (not bool(b["integrated_pass"]))})

    advantage_thresholds = [r["threshold"] for r in rows if r["uplift_advantage"]]
    out = {
        "packet": "P1181",
        "as_of": "2026-05-10",
        "thresholds": thresholds,
        "rows": rows,
        "uplift_advantage_thresholds": advantage_thresholds,
        "note": "Delta curve only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1181_multithreshold_delta_curve_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1181] thresholds={len(thresholds)} uplift_advantage={advantage_thresholds} wrote {out_path}")

if __name__ == "__main__":
    main()
