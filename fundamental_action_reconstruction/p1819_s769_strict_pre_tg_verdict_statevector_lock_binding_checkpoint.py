#!/usr/bin/env python3
"""P1819 S769 bind P1818 verdict/state-vector lock as mandatory pre-TG gate."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def main() -> None:
    p1818 = load("p1818_s768_strict_verdict_statevector_consistency_lock_checkpoint.json")
    p1813 = load("p1813_s763_strict_next_honest_step_bw_brst_cut_gate_execution_plan_checkpoint.json")
    p1814 = load("p1814_s764_strict_priority_closure_gap_normalization_checkpoint.json")
    p1815 = load("p1815_s765_strict_state_vector_priority_sync_checkpoint.json")

    lock_status = p1818.get("status", "OPEN_OBSTRUCTION_WITH_TRACE")
    lock_pass = lock_status == "LOCK_CONSISTENT"

    out = {
        "packet_id": "P1819",
        "stage_id": "S769",
        "status": "LOCK_CONSISTENT" if lock_pass else "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "binding": {
            "mandatory_pre_tg_gate": "P1818_S768",
            "required_status": "LOCK_CONSISTENT",
            "observed_status": lock_status,
            "enforced_scope": ["TG1_BW", "TG2_BRST", "TG3_CUT"],
        },
        "consistency_trace": {
            "p1813_active_gate": p1813.get("next_honest_step", {}).get("active_gate", "TG1_BW"),
            "p1814_tg_snapshot": p1814.get("theorem_gates", {}),
            "p1815_theorem_gates": p1815.get("state_vector", {}).get("theorem_gates", {}),
        },
        "decision": {
            "tg_gate_evaluation_allowed": lock_pass,
            "if_false_then": "KEEP_TG1_TG2_TG3_OPEN_WITH_TRACE",
            "pass_policy": "NO_FALSE_PASS_NO_AUTOPROMOTION",
        },
        "technical_progress": "Lock promoted from standalone check to mandatory theorem-gate precondition.",
        "proven": "TG evaluation chain now has explicit verdict/state-vector consistency guard dependency.",
        "open": "Global theorem gates remain OPEN pending real nonproxy witness residuals and QG checks.",
        "false_pass_risk": "Skipping this precheck could let stale SV artifacts leak into TG decisions.",
        "next_honest_step": "Execute real P1806 evidence runpack, rerun P1793->P1794->P1818, then evaluate TG1.",
        "lay_explanation": "Najpierw sprawdzamy, czy statusy są spójne; dopiero potem wolno oceniać główne bramki teorii.",
    }

    out_path = GEN / "p1819_s769_strict_pre_tg_verdict_statevector_lock_binding_checkpoint.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
