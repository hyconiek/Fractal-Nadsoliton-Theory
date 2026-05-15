#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"

# Recompute expected p1808 result directly from source checkpoints (same logic intent as emitter)
def detect(txt: str) -> str:
    if "W1 not FULL_EXPORT" in txt:
        return "NOT_FULL_EXPORT"
    if "W1 accepted as FULL_EXPORT" in txt:
        return "FULL_EXPORT_ACCEPTED"
    return "UNKNOWN"

src = {
    "p1779": detect((G/"p1779_s729_strict_current_priority_success_condition_tracker_checkpoint.json").read_text()),
    "p1782": detect((G/"p1782_s732_strict_priority_closure_gap_matrix_checkpoint.json").read_text()),
    "p1780": detect((G/"p1780_s730_strict_theorem_gate_freeze_until_joint_residual_witness_checkpoint.json").read_text()),
}
expected_conflict = len(set(src.values())) > 1
expected_tg1 = "OPEN_LOCKED_BY_W1_CONSISTENCY_CONFLICT" if expected_conflict else "OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL"

actual_p1808 = json.loads((G/"p1808_s758_strict_w1_full_export_gate_reconciliation_checkpoint.json").read_text())
actual_tg1 = actual_p1808.get("gate_vector",{}).get("TG1_BW")

out = {
  "packet_id":"P1811",
  "stage_id":"S761",
  "status":"PASS_ZERO_CONSISTENT" if actual_tg1 == expected_tg1 else "OPEN_OBSTRUCTION_WITH_TRACE",
  "expected_tg1_from_sources": expected_tg1,
  "actual_tg1_in_generated_p1808": actual_tg1,
  "source_states": src,
  "requires_regeneration": actual_tg1 != expected_tg1
}
(G/"p1811_s761_strict_generated_vs_emitter_consistency_audit_checkpoint.json").write_text(json.dumps(out, indent=2)+"\n")
print(str(G/"p1811_s761_strict_generated_vs_emitter_consistency_audit_checkpoint.json"))
