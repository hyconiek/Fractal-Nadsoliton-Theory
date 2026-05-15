#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1778 = GEN / "p1778_s728_strict_componentwise_h1_metric_joint_execution_lock_checkpoint.json"
IN1780 = GEN / "p1780_s730_strict_theorem_gate_freeze_until_joint_residual_witness_checkpoint.json"
IN1782 = GEN / "p1782_s732_strict_priority_closure_gap_matrix_checkpoint.json"
OUT = GEN / "p1783_s733_strict_priority_gate_exit_criteria_matrix_checkpoint.json"


def main() -> None:
    p1778 = json.loads(IN1778.read_text(encoding="utf-8"))
    p1780 = json.loads(IN1780.read_text(encoding="utf-8"))
    p1782 = json.loads(IN1782.read_text(encoding="utf-8"))

    joint_gate = p1778.get("gate_decision", {})
    freeze = p1780.get("theorem_gate_freeze", {})

    exit_matrix = {
        "W1_FULL_EXPORT": {
            "required": True,
            "met": joint_gate.get("joint_run_allowed_now", False),
            "source": "p1778.gate_decision",
        },
        "JOINT_H1_METRIC_COMPONENTWISE_RUN": {
            "required": True,
            "met": False,
            "source": "pending joint execution witness",
        },
        "JOINT_RESULT_CLASSIFIED_PASS_OR_OBSTRUCTION": {
            "required": True,
            "met": False,
            "allowed": ["PASS_ZERO", "OBSTRUCTION_WITH_DIVERGENCE_TRACE"],
        },
        "THEOREM_GATE_FREEZE_RELEASE": {
            "required": True,
            "met": not freeze.get("freeze_active", True),
            "source": "p1780.theorem_gate_freeze.freeze_active",
        },
    }

    payload = {
        "checkpoint": "P1783_S733",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1778_s728", "p1780_s730", "p1782_s732"],
        "priority_gate_exit_criteria_matrix": exit_matrix,
        "closure_gap_ref": p1782.get("priority_closure_gap_matrix", {}),
        "global_exit_ready": all(v.get("met", False) for v in exit_matrix.values()),
        "proofs_and_nonproofs": {
            "proven_now": [
                "A formal exit-criteria matrix for promoting priority gates is exported with machine-readable required/met flags.",
            ],
            "still_open": [
                "Exit remains not ready because joint witness and freeze-release criteria are unmet.",
            ],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Spełnić kolejno wszystkie kryteria exit-matrix (W1 FULL_EXPORT, joint run, klasyfikacja PASS/OBSTRUCTION, release freeze) i dopiero wtedy promować theorem gates.",
        "lay_summary": "To lista warunków wyjścia: dopóki choć jeden warunek nie jest spełniony, system uczciwie nie pozwala przejść do końcowych bramek.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
