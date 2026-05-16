#!/usr/bin/env python3
"""P1828 S778 strict TG1 blocker-closure dashboard checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def main() -> None:
    p1823 = load("p1823_s773_strict_tg1_preflight_gate_decision_checkpoint.json")
    p1826 = load("p1826_s776_strict_s1_to_tg1_readiness_projection_checkpoint.json")
    p1827 = load("p1827_s777_strict_s1_target_execution_order_checkpoint.json")

    preflight = p1823.get("tg1_preflight", {})
    progress = p1826.get("s1_delivery_progress", {})
    queue = p1827.get("ordered_s1_execution_targets", [])

    pending = int(progress.get("pending_targets", 0))
    tg1_allowed = bool(preflight.get("tg1_run_allowed", False))

    blocker_class = "NONE" if tg1_allowed else ("S1_ARTIFACT_DEFICIT" if pending > 0 else "PRECONDITION_MISMATCH")

    out = {
        "packet_id": "P1828",
        "stage_id": "S778",
        "status": "PASS_ZERO" if tg1_allowed else "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "tg1_dashboard": {
            "tg1_allowed": tg1_allowed,
            "blocker_class": blocker_class,
            "pending_targets": pending,
            "readiness_ratio": progress.get("readiness_ratio", 0.0),
            "lock_precondition_observed": preflight.get("lock_precondition_observed", "OPEN_UNKNOWN"),
            "s1_readiness_observed": preflight.get("s1_readiness_observed", "OPEN_UNKNOWN"),
        },
        "priority_execution_head": queue[:3],
        "technical_progress": "TG1 blockers are now classified into a single closure dashboard with explicit blocker class and execution head.",
        "proven": "Current blocker is classifiable as S1 artifact deficit under strict preflight logic.",
        "open": "TG1 remains blocked while pending target count is nonzero.",
        "false_pass_risk": "Without blocker-class tracking, teams may misreport readiness despite unresolved S1 deficits.",
        "next_honest_step": "Resolve queue head targets from P1827, update P1825 delivery tracker, then rerun P1826 and P1823.",
        "lay_explanation": "Mamy teraz jedną tablicę: pokazuje dokładnie, dlaczego TG1 stoi i które 3 zadania robić najpierw.",
    }

    path = GEN / "p1828_s778_strict_tg1_blocker_closure_dashboard_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
