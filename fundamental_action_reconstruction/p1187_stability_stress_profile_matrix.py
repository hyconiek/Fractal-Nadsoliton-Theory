#!/usr/bin/env python3
"""P1187 stability stress profile matrix across candidates and thresholds."""
from __future__ import annotations
import json, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def run_once(candidate: Path, threshold: float) -> bool:
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as tf:
        reg = Path(tf.name)
        tf.write(json.dumps({"candidates": [str(candidate)]}))
    cmd = [sys.executable, str(ROOT / "p1169_e2e_candidate_integration.py"), "--strict-e2e", "--require-out-of-locality-robustness", "--robustness-threshold", str(threshold), "--candidate", str(candidate), "--registry", str(reg)]
    subprocess.run(cmd, capture_output=True, text=True)
    s = json.loads((GEN / "p1169_e2e_candidate_integration_summary.json").read_text(encoding="utf-8"))
    return bool(s.get("integrated_pass"))


def main() -> None:
    candidates = {
        "baseline_neighbor_1": (GEN / "p1170_neighbor_1.json").resolve(),
        "uplift_neighbor_3": (GEN / "p1170_neighbor_3.json").resolve(),
    }
    thresholds = [0.55, 0.60, 0.65, 0.70]
    repeats = 3

    matrix = {}
    for cname, cpath in candidates.items():
        matrix[cname] = {}
        for t in thresholds:
            outcomes = [run_once(cpath, t) for _ in range(repeats)]
            matrix[cname][str(t)] = {
                "repeats": repeats,
                "pass_count": sum(1 for x in outcomes if x),
                "stability_rate": sum(1 for x in outcomes if x) / repeats,
                "outcomes": outcomes,
            }

    out = {
        "packet": "P1187",
        "as_of": "2026-05-10",
        "thresholds": thresholds,
        "repeats": repeats,
        "matrix": matrix,
        "note": "Stress profile only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1187_stability_stress_profile_matrix_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1187] wrote {out_path}")

if __name__ == "__main__":
    main()
