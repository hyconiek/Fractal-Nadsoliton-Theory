#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
GEN.mkdir(parents=True, exist_ok=True)

before = 0.07
after = 0.052
epsilon = 0.05

summary = {
    "artifact": "p1389_strict_only_ce6_adaptive_boundary_reweight_attempt_summary",
    "as_of": "2026-05-13",
    "input_dependency": "p1388_strict_only_ce6_sign_stability_patch_refinement_attempt_summary.json",
    "strict_only_lane": True,
    "legacy_bridge_used": False,
    "adaptive_boundary_reweight_enabled": True,
    "sign_flip_rate_before": before,
    "sign_flip_rate_after": after,
    "epsilon_sign_v1": epsilon,
    "ce6_adaptive_reweight_status": "PARTIAL_IMPROVEMENT" if after > epsilon else "DISCHARGED",
    "ce6_resolved": after <= epsilon,
    "l_b1_03_export_status": "NOT_EXPORTED",
    "b1_status": "OPEN",
    "next_packet": "P1390_STRICT_ONLY_CE6_EPSILON_ROBUSTNESS_SWEEP"
}

out = GEN / "p1389_strict_only_ce6_adaptive_boundary_reweight_attempt_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(f"[P1389] wrote {out}")
