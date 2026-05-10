#!/usr/bin/env python3
"""P1183 automated promotion gate based on strict artifacts."""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    p1180 = json.loads((GEN / "p1180_e2e_candidate_delta_comparison_summary.json").read_text(encoding="utf-8"))
    p1182 = json.loads((GEN / "p1182_robustness_band_lock_check_summary.json").read_text(encoding="utf-8"))
    p1152 = json.loads((GEN / "p1152_strict_candidate_registry_runner_summary.json").read_text(encoding="utf-8"))

    winner = (p1180.get("winner") or {}).get("candidate")
    band_lock_pass = bool(p1182.get("band_lock_pass"))
    shortlist_consistency_ok = bool(p1152.get("shortlist_consistency_enforced")) and int(p1152.get("shortlist_consistency_returncode", 1)) == 0
    safe_region_mode = bool(p1152.get("require_safe_region"))
    registry_pass = int(p1152.get("failed", 1)) == 0 and int(p1152.get("passed", 0)) >= 1

    promote = bool(winner) and band_lock_pass and shortlist_consistency_ok and safe_region_mode and registry_pass

    out = {
        "packet": "P1183",
        "as_of": "2026-05-10",
        "winner_candidate": winner,
        "band_lock_pass": band_lock_pass,
        "shortlist_consistency_ok": shortlist_consistency_ok,
        "safe_region_mode": safe_region_mode,
        "registry_pass": registry_pass,
        "promote_candidate": promote,
        "note": "Promotion gate only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1183_candidate_promotion_gate_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1183] promote_candidate={promote} winner={winner} wrote {out_path}")

if __name__ == "__main__":
    main()
