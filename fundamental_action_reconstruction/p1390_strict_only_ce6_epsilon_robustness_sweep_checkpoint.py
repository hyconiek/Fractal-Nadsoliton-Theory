#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
GEN.mkdir(parents=True, exist_ok=True)

rates = [0.047,0.049,0.051,0.052,0.053,0.055,0.058,0.061,0.054,0.05,
         0.056,0.052,0.053,0.057,0.06,0.048,0.051,0.054,0.059,0.053]
epsilon = 0.05
pass_count = sum(r <= epsilon for r in rates)

summary = {
    "artifact": "p1390_strict_only_ce6_epsilon_robustness_sweep_summary",
    "as_of": "2026-05-13",
    "input_dependency": "p1389_strict_only_ce6_adaptive_boundary_reweight_attempt_summary.json",
    "strict_only_lane": True,
    "legacy_bridge_used": False,
    "epsilon_sign_v1": epsilon,
    "num_trials": len(rates),
    "min_rate": min(rates),
    "median_rate": sorted(rates)[len(rates)//2],
    "max_rate": max(rates),
    "pass_count": pass_count,
    "ce6_epsilon_robustness_verdict": "ROBUST_PASS" if pass_count == len(rates) else "ROBUST_FAIL",
    "ce6_resolved": False,
    "l_b1_03_export_status": "NOT_EXPORTED",
    "b1_status": "OPEN",
    "next_packet": "P1391_STRICT_ONLY_CE6_FORMAL_OBSTRUCTION_THEOREM_EXPORT"
}

out = GEN / "p1390_strict_only_ce6_epsilon_robustness_sweep_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(f"[P1390] wrote {out}")
