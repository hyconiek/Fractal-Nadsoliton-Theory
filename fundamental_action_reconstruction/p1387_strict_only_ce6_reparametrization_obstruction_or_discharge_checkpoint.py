#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
GEN.mkdir(parents=True, exist_ok=True)

summary = {
    "artifact": "p1387_strict_only_ce6_reparametrization_obstruction_or_discharge_summary",
    "as_of": "2026-05-13",
    "input_dependency": "p1386_strict_only_counterexample_ce4_ce6_resolution_attempt_summary.json",
    "strict_only_lane": True,
    "legacy_bridge_used": False,
    "target_counterexample": "ce6",
    "ce6_tests_v1": {
        "overlap_reparam_perturbation": "performed",
        "selector_drift_check": "performed",
        "boundary_sign_stability": "failed_subregion_detected"
    },
    "ce6_verdict": "OBSTRUCTION_EXPORTED_V1",
    "theorem_discharge": "NOT_ACHIEVED",
    "l_b1_03_export_status": "NOT_EXPORTED",
    "b1_status": "OPEN",
    "next_packet": "P1388_STRICT_ONLY_CE6_SIGN_STABILITY_PATCH_REFINEMENT_ATTEMPT"
}

out = GEN / "p1387_strict_only_ce6_reparametrization_obstruction_or_discharge_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(f"[P1387] wrote {out}")
