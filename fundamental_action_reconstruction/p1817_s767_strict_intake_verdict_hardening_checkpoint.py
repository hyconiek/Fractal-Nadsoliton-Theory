#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

p1793 = json.loads((GEN / "p1793_s743_sv1_sv5_evidence_intake_validation_result.json").read_text())
p1816 = json.loads((GEN / "p1816_s766_strict_unified_runpack_intake_validation_sync_checkpoint.json").read_text())

out = {
    "packet_id": "P1817",
    "stage_id": "S767",
    "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    "hardening_applied": "P1793_PASS_ZERO_requires_all_ledgers_PASS_ZERO",
    "post_fix_validation": {
        "p1793_verdict": p1793.get("verdict", "OPEN_OBSTRUCTION_WITH_TRACE"),
        "p1793_obstruction_tags": p1793.get("obstruction_tags", [])
    },
    "state_vector_policy": p1816.get("gate_policy", {}),
    "false_pass_guard": "SCHEMA_VALIDATION_CANNOT_PROMOTE_THEOREM_GATES",
    "next_honest_step": "inject_real_P1806_execution_ledger_and_rerun_validation_chain"
}

path = GEN / "p1817_s767_strict_intake_verdict_hardening_checkpoint.json"
path.write_text(json.dumps(out, indent=2) + "\n")
print(path)
