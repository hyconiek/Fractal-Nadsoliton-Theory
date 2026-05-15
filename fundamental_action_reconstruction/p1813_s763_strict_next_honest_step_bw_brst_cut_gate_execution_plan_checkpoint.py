#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"
P1812 = json.loads((G/"p1812_s762_strict_canonical_gate_status_source_checkpoint.json").read_text())
P1806 = json.loads((G/"p1806_s756_strict_unified_nonproxy_ea_eh_elg_residual_runpack_contract_checkpoint.json").read_text())

tg1 = P1812.get("canonical_gate_status", {}).get("TG1_BW", "OPEN_UNKNOWN")
runpack_ready = bool(P1806.get("status", "").startswith("PASS_ZERO"))

if tg1 == "PASS_ZERO":
    priority = "ADVANCE_TO_TG2_BRST_NILPOTENCY"
    next_gate = "TG2_BRST"
else:
    priority = "EXECUTE_UNIFIED_NONPROXY_BW_RESIDUAL_WITNESS"
    next_gate = "TG1_BW"

out = {
  "packet_id": "P1813",
  "stage_id": "S763",
  "status": "OPEN_OBSTRUCTION_WITH_TRACE" if tg1 != "PASS_ZERO" else "PASS_ZERO_READY_FOR_TG2",
  "inputs": {
    "p1812_tg1_bw": tg1,
    "p1806_runpack_contract_status": P1806.get("status", "UNKNOWN")
  },
  "classification": {
    "route": "strict_only",
    "scale": "global_nonproxy",
    "pass_policy": "binary_witness_only"
  },
  "next_honest_step": {
    "priority": priority,
    "active_gate": next_gate,
    "runpack_ready": runpack_ready,
    "required_output": ["PASS_ZERO", "OBSTRUCTION_WITH_TRACE"],
    "sync_order_after_verdict": ["P1810", "P1811", "P1812"]
  }
}

path = G/"p1813_s763_strict_next_honest_step_bw_brst_cut_gate_execution_plan_checkpoint.json"
path.write_text(json.dumps(out, indent=2) + "\n")
print(str(path))
