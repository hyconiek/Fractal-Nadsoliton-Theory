#!/usr/bin/env python3
"""P1185 post-certification drift check."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    return {"cmd": cmd, "returncode": p.returncode, "stdout": p.stdout.strip(), "stderr": p.stderr.strip()}


def main() -> None:
    reg = (GEN / "p1184_promoted_candidate_registry.json").resolve()
    r1152 = run([sys.executable, str(ROOT / "p1152_strict_candidate_registry_runner.py"), "--require-safe-region", "--rank-by-safe-margin", "--enforce-shortlist-consistency", str(reg)])
    s1152 = json.loads((GEN / "p1152_strict_candidate_registry_runner_summary.json").read_text(encoding="utf-8"))

    r1177 = run([sys.executable, str(ROOT / "p1177_shortlist_consistency_check.py"), str(GEN / "p1152_strict_candidate_registry_runner_summary.json")])
    s1177 = json.loads((GEN / "p1177_shortlist_consistency_check_summary.json").read_text(encoding="utf-8"))

    r1182 = run([sys.executable, str(ROOT / "p1182_robustness_band_lock_check.py")])
    s1182 = json.loads((GEN / "p1182_robustness_band_lock_check_summary.json").read_text(encoding="utf-8"))

    stable = (r1152["returncode"] == 0 and r1177["returncode"] == 0 and r1182["returncode"] == 0
              and int(s1152.get("failed", 1)) == 0
              and bool(s1177.get("overall_pass"))
              and bool(s1182.get("band_lock_pass")))

    out = {
        "packet": "P1185",
        "as_of": "2026-05-10",
        "registry": str(reg),
        "runs": {"p1152": r1152, "p1177": r1177, "p1182": r1182},
        "checks": {
            "p1152_failed": s1152.get("failed"),
            "p1177_overall_pass": s1177.get("overall_pass"),
            "p1182_band_lock_pass": s1182.get("band_lock_pass"),
        },
        "post_certification_stable": stable,
        "note": "Drift check only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1185_post_certification_drift_check_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1185] post_certification_stable={stable} wrote {out_path}")

if __name__ == "__main__":
    main()
