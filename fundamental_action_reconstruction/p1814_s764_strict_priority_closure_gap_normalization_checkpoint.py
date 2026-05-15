#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

def load(name: str) -> dict:
    return json.loads((GEN / name).read_text())

p1764 = load("p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json")
p1765 = load("p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json")
p1762 = load("p1762_s712_strict_boundary_control_contract_finalization_checkpoint.json")
p1763 = load("p1763_s713_strict_nonproxy_h1_4d_first_execution_attempt_checkpoint.json")
p1766 = load("p1766_s716_strict_forward_reverse_state_vector_update_with_bianchi_ward_gate_checkpoint.json")
p1801 = load("p1801_s751_strict_brst_nilpotency_witness_intake_and_gate_link_checkpoint.json")
p1802 = load("p1802_s752_strict_cutkosky_unitarity_witness_intake_and_tg3_gate_link_checkpoint.json")
p1812 = load("p1812_s762_strict_canonical_gate_status_source_checkpoint.json")

TG = p1812.get("canonical_gate_status", {})

def open_if_not_pass(value: str) -> str:
    return "PASS_ZERO" if value == "PASS_ZERO" else "OPEN_OBSTRUCTION_WITH_TRACE"

out = {
    "packet_id": "P1814",
    "stage_id": "S764",
    "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    "route": "strict_only",
    "priority_snapshot": {
        "explicit_nonproxy_E_A_mu": p1764.get("status", "OPEN"),
        "explicit_nonproxy_E_H": p1764.get("status", "OPEN"),
        "metric_EL_g_export": p1765.get("status", "OPEN"),
        "boundary_term_control": p1762.get("status", "OPEN"),
        "H1_4D_weak_form_readiness": p1763.get("status", "OPEN"),
        "bianchi_ward_consistency": p1766.get("status", "OPEN")
    },
    "theorem_gates": {
        "TG1_BW": open_if_not_pass(TG.get("TG1_BW", "OPEN")),
        "TG2_BRST": "OPEN_OBSTRUCTION_WITH_TRACE",
        "TG3_CUT": "OPEN_OBSTRUCTION_WITH_TRACE"
    },
    "evidence_requirements": {
        "TG1_required": "unified_nonproxy_BW_residual_witness_trace",
        "TG2_required": p1801.get("tg2_pass_requirements", []),
        "TG3_required": p1802.get("tg3_pass_requirements", [])
    },
    "no_false_pass_policy": {
        "pass_allowed_only_with": ["PASS_ZERO", "witness_trace", "residual_check"],
        "current_decision": "KEEP_ALL_REMAINING_GATES_OPEN"
    },
    "next_honest_step": "execute_P1806_unified_nonproxy_runpack_then_sync_P1810_P1811_P1812"
}

path = GEN / "p1814_s764_strict_priority_closure_gap_normalization_checkpoint.json"
path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
print(path)
