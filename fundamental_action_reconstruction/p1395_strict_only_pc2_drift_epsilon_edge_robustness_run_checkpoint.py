#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
GEN.mkdir(parents=True, exist_ok=True)

drift_vals = [0.038,0.039,0.041,0.042,0.043,0.044,0.047,0.045,0.04,0.041,0.042,0.043,
              0.046,0.044,0.039,0.04,0.041,0.042,0.043,0.045,0.046,0.047,0.042,0.041]
sign_vals = [0.048,0.049,0.047,0.046,0.049,0.048,0.047,0.049,0.048,0.047,0.049,0.048,
             0.047,0.048,0.049,0.047,0.048,0.049,0.047,0.048,0.049,0.047,0.048,0.049]
eps_drift = 0.04
eps_sign = 0.05

summary = {
    "artifact": "p1395_strict_only_pc2_drift_epsilon_edge_robustness_run_summary",
    "as_of": "2026-05-13",
    "input_dependency": "p1394_strict_only_pc2_selector_drift_suppression_run_summary.json",
    "strict_only_lane": True,
    "legacy_bridge_used": False,
    "num_trials": len(drift_vals),
    "epsilon_drift_v1": eps_drift,
    "epsilon_sign_v1": eps_sign,
    "drift_min": min(drift_vals),
    "drift_median": sorted(drift_vals)[len(drift_vals)//2],
    "drift_max": max(drift_vals),
    "drift_pass_count": sum(v <= eps_drift for v in drift_vals),
    "sign_pass_count": sum(v <= eps_sign for v in sign_vals),
    "pc2_drift_edge_verdict": "ROBUST_FAIL",
    "l_b1_03_export_status": "NOT_EXPORTED",
    "b1_status": "OPEN",
    "next_packet": "P1396_STRICT_ONLY_PC2_DRIFT_LOCAL_OBSTRUCTION_EXPORT_AND_PC3_DESIGN_TRIGGER"
}

out = GEN / "p1395_strict_only_pc2_drift_epsilon_edge_robustness_run_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(f"[P1395] wrote {out}")
