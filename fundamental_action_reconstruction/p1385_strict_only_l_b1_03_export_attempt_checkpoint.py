#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
GEN.mkdir(parents=True, exist_ok=True)

counterexamples = {
    "ce1": "resolved",
    "ce2": "resolved",
    "ce3": "resolved",
    "ce4": "unresolved",
    "ce5": "resolved",
    "ce6": "unresolved",
}
resolved = sum(v == "resolved" for v in counterexamples.values())
total = len(counterexamples)

summary = {
    "artifact": "p1385_strict_only_l_b1_03_export_attempt_summary",
    "as_of": "2026-05-13",
    "input_dependency": "p1384_strict_only_selector_compatibility_discharge_attempt_summary.json",
    "strict_only_lane": True,
    "legacy_bridge_used": False,
    "target_lemma": "L-B1-03",
    "counterexample_set_v1": counterexamples,
    "resolved_counterexamples": resolved,
    "total_counterexamples": total,
    "global_selector_compatibility_bound_verified": False,
    "l_b1_03_export_status": "PARTIAL_NOT_EXPORTED",
    "g3_status": "PARTIAL_QW2191_NARROWED",
    "b1_status": "OPEN",
    "next_packet": "P1386_STRICT_ONLY_COUNTEREXAMPLE_CE4_CE6_RESOLUTION_ATTEMPT"
}

out = GEN / "p1385_strict_only_l_b1_03_export_attempt_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(f"[P1385] wrote {out}")
