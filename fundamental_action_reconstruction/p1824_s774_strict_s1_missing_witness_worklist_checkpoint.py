#!/usr/bin/env python3
"""P1824 S774 strict S1 missing-witness worklist checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def workitem(target: str) -> dict:
    return {
        "target": target,
        "required_artifacts": [
            f"{target}::explicit_covariant_formula_export",
            f"{target}::residual_expression_export",
            f"{target}::witness_trace_log",
            f"{target}::independent_replay_check",
        ],
        "acceptance": "all artifacts present and status upgraded to PASS_ZERO by source checkpoint",
    }


def main() -> None:
    p1822 = load("p1822_s772_strict_s1_evidence_pack_readiness_checkpoint.json")
    p1823 = load("p1823_s773_strict_tg1_preflight_gate_decision_checkpoint.json")

    missing = p1822.get("missing_witness_targets", [])
    tg1_allowed = p1823.get("tg1_preflight", {}).get("tg1_run_allowed", False)

    worklist = [workitem(t) for t in missing]

    out = {
        "packet_id": "P1824",
        "stage_id": "S774",
        "status": "PASS_ZERO" if not missing else "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "input_gate": {
            "p1822_status": p1822.get("status", "OPEN_UNKNOWN"),
            "p1823_tg1_run_allowed": tg1_allowed,
        },
        "missing_count": len(missing),
        "s1_missing_witness_worklist": worklist,
        "execution_rule": "complete every workitem before rerunning P1822/P1823",
        "technical_progress": "Missing S1 evidence obligations translated into explicit per-target artifact worklist.",
        "proven": "TG1 remains blocked until each missing target has formula+residual+witness+replay artifacts.",
        "open": "Worklist is non-empty, so strict TG1 preflight cannot pass yet.",
        "false_pass_risk": "Claiming readiness without per-target replayable artifacts would create non-auditable PASS drift.",
        "next_honest_step": "Produce artifacts for each workitem target, then rerun P1822 and P1823.",
        "lay_explanation": "Rozbijamy brakujące dowody na konkretne zadania; dopóki lista nie jest pusta, TG1 nie startuje.",
    }

    path = GEN / "p1824_s774_strict_s1_missing_witness_worklist_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
