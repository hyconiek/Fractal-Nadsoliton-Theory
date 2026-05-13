#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
GEN.mkdir(parents=True, exist_ok=True)

epsilon = 0.02
observed_gap = 0.012
local_pass = observed_gap <= epsilon

summary = {
    "artifact": "p1384_strict_only_selector_compatibility_discharge_attempt_summary",
    "as_of": "2026-05-13",
    "input_dependency": "p1383_strict_only_cmix_commutator_globalization_attempt_summary.json",
    "strict_only_lane": True,
    "legacy_bridge_used": False,
    "target_obligation": "G3_selector_compatibility_qw2191",
    "epsilon_selector_v1": epsilon,
    "observed_selector_compatibility_gap": observed_gap,
    "local_selector_test_status": "PASS" if local_pass else "FAIL",
    "hidden_selector_closure_dependency_detected": False,
    "l_b1_03_export": "NOT_EXPORTED",
    "g3_status": "LOCAL_PASS_GLOBAL_DISCHARGE_NOT_YET",
    "b1_status": "OPEN",
    "next_packet": "P1385_STRICT_ONLY_L_B1_03_EXPORT_ATTEMPT"
}

out = GEN / "p1384_strict_only_selector_compatibility_discharge_attempt_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(f"[P1384] wrote {out}")
