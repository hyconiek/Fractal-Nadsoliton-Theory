#!/usr/bin/env python3
"""P1830 S780 strict TG1 governance gate checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def main() -> None:
    p1828 = load("p1828_s778_strict_tg1_blocker_closure_dashboard_checkpoint.json")
    p1829 = load("p1829_s779_strict_tg1_blocker_burndown_checkpoint.json")

    dash = p1828.get("tg1_dashboard", {})
    burn = p1829.get("burndown", {})

    blocker_class = dash.get("blocker_class", "UNKNOWN")
    pending = int(burn.get("pending_targets", 0))
    burndown_ratio = float(burn.get("burndown_ratio", 0.0))
    tg1_allowed = bool(dash.get("tg1_allowed", False))

    governance_state = "BLOCKED_BY_S1_EVIDENCE" if blocker_class == "S1_ARTIFACT_DEFICIT" else "REVIEW_REQUIRED"
    if tg1_allowed and pending == 0:
        governance_state = "READY_FOR_TG1_EXECUTION"

    out = {
        "packet_id": "P1830",
        "stage_id": "S780",
        "status": "PASS_ZERO" if governance_state == "READY_FOR_TG1_EXECUTION" else "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "governance_gate": {
            "state": governance_state,
            "blocker_class": blocker_class,
            "pending_targets": pending,
            "burndown_ratio": burndown_ratio,
            "tg1_allowed_flag": tg1_allowed,
        },
        "release_rule": "TG1 execution release only if READY_FOR_TG1_EXECUTION",
        "technical_progress": "TG1 blocker status is now translated into a governance gate state for release/no-release decisions.",
        "proven": "Current governance state remains blocked while S1 evidence deficit persists.",
        "open": "No release permission for TG1 execution under current blocker state.",
        "false_pass_risk": "Executing TG1 outside governance gate can bypass strict evidence closure and create false theorem momentum.",
        "next_honest_step": "Close S1 deficits until governance state flips to READY_FOR_TG1_EXECUTION, then run TG1 once.",
        "lay_explanation": "To formalna zgoda/odmowa startu TG1: dziś status mówi 'blokada', więc nie uruchamiamy bramki.",
    }

    path = GEN / "p1830_s780_strict_tg1_governance_gate_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
