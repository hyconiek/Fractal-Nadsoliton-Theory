#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

def load(n: str) -> dict:
    return json.loads((GEN / n).read_text())

p1812 = load("p1812_s762_strict_canonical_gate_status_source_checkpoint.json")
p1814 = load("p1814_s764_strict_priority_closure_gap_normalization_checkpoint.json")
p1801 = load("p1801_s751_strict_brst_nilpotency_witness_intake_and_gate_link_checkpoint.json")
p1802 = load("p1802_s752_strict_cutkosky_unitarity_witness_intake_and_tg3_gate_link_checkpoint.json")

canonical = p1812.get("canonical_gate_status", {})

tg1 = canonical.get("TG1_BW", "OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL")
tg2 = canonical.get("TG2_BRST", "LOCKED_BY_TG1")
tg3 = canonical.get("TG3_CUT", "LOCKED_BY_TG2")

state = {
    "forward_chain": "EXPORTED_NONPROXY_KEEP_OPEN",
    "reverse_chain": {
        "H1_4D": p1814.get("priority_snapshot", {}).get("H1_4D_weak_form_readiness", "OPEN"),
        "metric_residual": p1814.get("theorem_gates", {}).get("TG1_BW", "OPEN_OBSTRUCTION_WITH_TRACE"),
        "bianchi_ward": p1814.get("priority_snapshot", {}).get("bianchi_ward_consistency", "OPEN")
    },
    "theorem_gates": {
        "TG1_BW": tg1,
        "TG2_BRST": tg2,
        "TG3_CUT": tg3
    }
}

out = {
    "packet_id": "P1815",
    "stage_id": "S765",
    "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    "route": "strict_only",
    "state_vector": state,
    "hard_requirements": {
        "TG2_BRST": p1801.get("tg2_pass_requirements", []),
        "TG3_CUT": p1802.get("tg3_pass_requirements", [])
    },
    "false_pass_guard": "NO_PROMOTION_WITHOUT_TG1_BINARY_WITNESS",
    "next_honest_step": "run_P1806_and_publish_binary_TG1_verdict_with_trace"
}

path = GEN / "p1815_s765_strict_state_vector_priority_sync_checkpoint.json"
path.write_text(json.dumps(out, indent=2) + "\n")
print(path)
