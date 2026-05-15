#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

def load(n: str) -> dict:
    return json.loads((GEN / n).read_text())

p1815 = load("p1815_s765_strict_state_vector_priority_sync_checkpoint.json")
validation = load("p1816_s766_unified_runpack_intake_validation_result.json")

val_verdict = validation.get("verdict", "OPEN_OBSTRUCTION_WITH_TRACE")

out = {
    "packet_id": "P1816",
    "stage_id": "S766",
    "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    "route": "strict_only",
    "input_validation": {
        "validator": "P1793",
        "validation_verdict": val_verdict,
        "obstruction_tags": validation.get("obstruction_tags", []),
        "ledger_count": validation.get("ledger_count", 0)
    },
    "state_vector_sync": p1815.get("state_vector", {}),
    "gate_policy": {
        "TG1_BW": "NO_PROMOTION_WITHOUT_GLOBAL_NONPROXY_RESIDUAL_WITNESS",
        "TG2_BRST": "LOCKED_BY_TG1",
        "TG3_CUT": "LOCKED_BY_TG2"
    },
    "false_pass_guard": "INPUT_VALIDATION_IS_NOT_PHYSICS_CLOSURE",
    "next_honest_step": "replace_demo_intake_with_real_P1806_evidence_and_rerun_P1793_then_P1794"
}

path = GEN / "p1816_s766_strict_unified_runpack_intake_validation_sync_checkpoint.json"
path.write_text(json.dumps(out, indent=2) + "\n")
print(path)
