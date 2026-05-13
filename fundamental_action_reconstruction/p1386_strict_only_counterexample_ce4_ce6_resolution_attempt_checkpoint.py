#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
GEN.mkdir(parents=True, exist_ok=True)

ce_status = {"ce4": "resolved", "ce6": "unresolved"}
all_resolved = all(v == "resolved" for v in ce_status.values())

summary = {
    "artifact": "p1386_strict_only_counterexample_ce4_ce6_resolution_attempt_summary",
    "as_of": "2026-05-13",
    "input_dependency": "p1385_strict_only_l_b1_03_export_attempt_summary.json",
    "strict_only_lane": True,
    "legacy_bridge_used": False,
    "target_counterexamples": ["ce4", "ce6"],
    "counterexample_status": ce_status,
    "ce_resolution_status": "FULL" if all_resolved else "PARTIAL",
    "l_b1_03_export_status": "EXPORTED" if all_resolved else "NOT_EXPORTED",
    "g3_status": "QW2191_DISCHARGED" if all_resolved else "PARTIAL_QW2191_NARROWED",
    "b1_status": "OPEN",
    "next_packet": "P1387_STRICT_ONLY_CE6_REPARAMETRIZATION_OBSTRUCTION_OR_DISCHARGE"
}

out = GEN / "p1386_strict_only_counterexample_ce4_ce6_resolution_attempt_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(f"[P1386] wrote {out}")
