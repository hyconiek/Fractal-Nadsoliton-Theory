#!/usr/bin/env python3
"""P1829 S779 strict TG1 blocker burndown checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def main() -> None:
    p1828 = load("p1828_s778_strict_tg1_blocker_closure_dashboard_checkpoint.json")
    p1825 = load("p1825_s775_strict_s1_worklist_delivery_tracker_checkpoint.json")

    tg1 = p1828.get("tg1_dashboard", {})
    summary = p1825.get("summary", {})

    total = int(summary.get("total_targets", 0))
    pending = int(summary.get("pending_targets", total))
    ready = int(summary.get("ready_targets", 0))

    blocker_burndown_ratio = 0.0 if total == 0 else (total - pending) / total

    out = {
        "packet_id": "P1829",
        "stage_id": "S779",
        "status": "PASS_ZERO" if pending == 0 and tg1.get("tg1_allowed", False) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "burndown": {
            "blocker_class": tg1.get("blocker_class", "UNKNOWN"),
            "total_targets": total,
            "ready_targets": ready,
            "pending_targets": pending,
            "burndown_ratio": blocker_burndown_ratio
        },
        "gate_projection": {
            "tg1_allowed_now": tg1.get("tg1_allowed", False),
            "condition_for_go": "pending_targets == 0 and preflight == PASS_ZERO"
        },
        "technical_progress": "Blocker dashboard is extended with a quantitative burndown metric for TG1 closure tracking.",
        "proven": "As long as pending_targets > 0, TG1 go condition remains unsatisfied.",
        "open": "S1 blocker burndown is incomplete.",
        "false_pass_risk": "Narrative progress without burndown closure can mask unresolved critical blockers.",
        "next_honest_step": "Burn down pending targets to zero via P1825 updates, then rerun P1823 and P1828.",
        "lay_explanation": "To licznik spalania blokad: dopóki nie zejdzie do zera, TG1 nie przechodzi.",
    }

    path = GEN / "p1829_s779_strict_tg1_blocker_burndown_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
