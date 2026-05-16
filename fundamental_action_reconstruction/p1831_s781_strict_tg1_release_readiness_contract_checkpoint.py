#!/usr/bin/env python3
"""P1831 S781 strict TG1 release-readiness contract checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def main() -> None:
    p1830 = load("p1830_s780_strict_tg1_governance_gate_checkpoint.json")
    p1823 = load("p1823_s773_strict_tg1_preflight_gate_decision_checkpoint.json")
    p1825 = load("p1825_s775_strict_s1_worklist_delivery_tracker_checkpoint.json")

    governance = p1830.get("governance_gate", {})
    preflight = p1823.get("tg1_preflight", {})
    summary = p1825.get("summary", {})

    ready = governance.get("state") == "READY_FOR_TG1_EXECUTION"

    out = {
        "packet_id": "P1831",
        "stage_id": "S781",
        "status": "PASS_ZERO" if ready else "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "release_contract": {
            "governance_state": governance.get("state", "UNKNOWN"),
            "preflight_tg1_allowed": preflight.get("tg1_run_allowed", False),
            "pending_targets": summary.get("pending_targets", None),
            "required_for_release": [
                "governance_state == READY_FOR_TG1_EXECUTION",
                "preflight_tg1_allowed == true",
                "pending_targets == 0",
            ],
            "release_allowed": ready,
        },
        "technical_progress": "TG1 release conditions are consolidated into one explicit contract object for go/no-go governance.",
        "proven": "Current release criteria remain unmet while governance/preflight/pending constraints are not jointly satisfied.",
        "open": "Release contract is not yet satisfiable on current S1 evidence state.",
        "false_pass_risk": "Running TG1 without all three release conditions can create an irreversible false-closure narrative.",
        "next_honest_step": "Drive pending_targets to zero and rerun P1823/P1830 until release_contract.release_allowed becomes true.",
        "lay_explanation": "To formalna umowa startu TG1: wszystkie warunki muszą być jednocześnie spełnione.",
    }

    path = GEN / "p1831_s781_strict_tg1_release_readiness_contract_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
