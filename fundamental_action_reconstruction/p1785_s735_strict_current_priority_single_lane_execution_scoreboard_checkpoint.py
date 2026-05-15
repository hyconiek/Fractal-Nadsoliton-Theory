#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1783 = GEN / "p1783_s733_strict_priority_gate_exit_criteria_matrix_checkpoint.json"
IN1784 = GEN / "p1784_s734_strict_current_priority_blocker_to_action_contract_checkpoint.json"
OUT = GEN / "p1785_s735_strict_current_priority_single_lane_execution_scoreboard_checkpoint.json"


def main() -> None:
    p1783 = json.loads(IN1783.read_text(encoding="utf-8"))
    p1784 = json.loads(IN1784.read_text(encoding="utf-8"))

    exit_matrix = p1783.get("priority_gate_exit_criteria_matrix", {})

    scoreboard = [
        {
            "lane_step": "L1_A1_W1_FULL_EXPORT",
            "status": "OPEN",
            "ready_to_execute": True,
            "completion_condition": "W1_FULL_EXPORT met=True",
            "source": exit_matrix.get("W1_FULL_EXPORT", {}),
        },
        {
            "lane_step": "L2_A2_JOINT_H1_METRIC_RUN",
            "status": "OPEN_BLOCKED_BY_L1",
            "ready_to_execute": False,
            "completion_condition": "JOINT_H1_METRIC_COMPONENTWISE_RUN met=True + result classified PASS_ZERO/OBSTRUCTION_WITH_DIVERGENCE_TRACE",
            "source": exit_matrix.get("JOINT_H1_METRIC_COMPONENTWISE_RUN", {}),
        },
        {
            "lane_step": "L3_A3_BW_BRST_CUT_PROGRESS",
            "status": "OPEN_BLOCKED_BY_L1_L2",
            "ready_to_execute": False,
            "completion_condition": "THEOREM_GATE_FREEZE_RELEASE met=True",
            "source": exit_matrix.get("THEOREM_GATE_FREEZE_RELEASE", {}),
        },
    ]

    payload = {
        "checkpoint": "P1785_S735",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1783_s733", "p1784_s734"],
        "single_lane_execution_scoreboard": scoreboard,
        "global_exit_ready": p1783.get("global_exit_ready", False),
        "proofs_and_nonproofs": {
            "proven_now": [
                "Single-lane ordering for current-priority execution is explicitly published with step-level readiness flags.",
            ],
            "still_open": [
                "Only L1 is executable now; L2 and L3 remain blocked by unmet upstream conditions.",
            ],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Wykonać L1 (A1), po sukcesie odblokować L2, a następnie L3 — bez zmiany kolejności i bez promotion theorem-gates przed witnessami.",
        "lay_summary": "To tablica etapów: pokazuje, co można zrobić teraz, a co musi poczekać na wcześniejsze kroki.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
