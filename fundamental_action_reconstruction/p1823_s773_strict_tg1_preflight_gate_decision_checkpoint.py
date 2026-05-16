#!/usr/bin/env python3
"""P1823 S773 strict TG1 preflight gate decision checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def main() -> None:
    p1819 = load("p1819_s769_strict_pre_tg_verdict_statevector_lock_binding_checkpoint.json")
    p1821 = load("p1821_s771_strict_priority_lane_sequencer_checkpoint.json")
    p1822 = load("p1822_s772_strict_s1_evidence_pack_readiness_checkpoint.json")

    lock_ok = p1819.get("status") == "LOCK_CONSISTENT"
    s1_ready = p1822.get("status") == "PASS_ZERO"

    lane = p1821.get("execution_lane", [])
    s2 = next((x for x in lane if x.get("step") == "S2_EXECUTE_TG1_UNIFIED_NONPROXY_RESIDUAL_RUN"), {})

    can_run_tg1 = lock_ok and s1_ready

    out = {
        "packet_id": "P1823",
        "stage_id": "S773",
        "status": "PASS_ZERO" if can_run_tg1 else "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "tg1_preflight": {
            "lock_precondition_required": "P1819 == LOCK_CONSISTENT",
            "lock_precondition_observed": p1819.get("status", "OPEN_UNKNOWN"),
            "s1_readiness_required": "P1822 == PASS_ZERO",
            "s1_readiness_observed": p1822.get("status", "OPEN_UNKNOWN"),
            "s1_missing_witness_targets": p1822.get("missing_witness_targets", []),
            "tg1_run_allowed": can_run_tg1,
        },
        "s2_contract": {
            "step": s2.get("step", "S2_EXECUTE_TG1_UNIFIED_NONPROXY_RESIDUAL_RUN"),
            "required_output": s2.get("required_output", []),
            "completion_rule": s2.get("completion_rule", "TG1 verdict exported with witness trace"),
        },
        "technical_progress": "TG1 entry is now guarded by a two-key preflight gate: lock consistency + S1 evidence readiness.",
        "proven": "TG1 cannot be launched under strict policy when S1 evidence pack is incomplete.",
        "open": "S1 witness targets remain incomplete, so TG1 run stays blocked.",
        "false_pass_risk": "Running TG1 without this preflight could create ambiguous non-reproducible partial-pass narratives.",
        "next_honest_step": "Complete all P1822 missing witness targets, rerun P1822 then P1823, and only then execute TG1 runpack.",
        "lay_explanation": "To jest zielone światło przed TG1: potrzebujesz i spójnego locka, i kompletnej teczki dowodów S1.",
    }

    path = GEN / "p1823_s773_strict_tg1_preflight_gate_decision_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
