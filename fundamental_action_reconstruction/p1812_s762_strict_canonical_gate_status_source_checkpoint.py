#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"
P1810 = json.loads((G/"p1810_s760_strict_tg1_lock_precedence_resolution_checkpoint.json").read_text())
P1811 = json.loads((G/"p1811_s761_strict_generated_vs_emitter_consistency_audit_checkpoint.json").read_text())

regeneration_required = bool(P1811.get("requires_regeneration", False))
tg1 = P1810.get("tg1_resolved", "OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL")

if regeneration_required:
    status = "OPEN_OBSTRUCTION_WITH_TRACE"
    canonical_tg1 = "OPEN_LOCKED_BY_ARTIFACT_REGENERATION_REQUIRED"
else:
    canonical_tg1 = tg1
    if canonical_tg1 == "PASS_ZERO":
        status = "PASS_ZERO_CANONICAL_SOURCE_READY"
    else:
        status = "OPEN_OBSTRUCTION_WITH_TRACE"

tg2 = "OPEN_REQUIRES_NILPOTENCY_WITNESS" if canonical_tg1 == "PASS_ZERO" else "LOCKED_BY_TG1"
tg3 = "OPEN_REQUIRES_UNITARITY_WITNESS" if tg2 == "PASS_ZERO" else "LOCKED_BY_TG2"

out = {
  "packet_id":"P1812",
  "stage_id":"S762",
  "status": status,
  "inputs": {
    "p1810_tg1_resolved": tg1,
    "p1811_requires_regeneration": regeneration_required
  },
  "canonical_gate_status": {
    "TG1_BW": canonical_tg1,
    "TG2_BRST": tg2,
    "TG3_CUT": tg3
  },
  "usage_policy": "USE_ONLY_THIS_ARTIFACT_FOR_DOWNSTREAM_GATE_STATE"
}

path = G/"p1812_s762_strict_canonical_gate_status_source_checkpoint.json"
path.write_text(json.dumps(out, indent=2)+"\n")
print(str(path))
