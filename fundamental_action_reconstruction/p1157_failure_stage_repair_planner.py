#!/usr/bin/env python3
"""P1157 failure-stage repair planner.

Transforms ranking with failure diagnostics into actionable repair templates.
Methodological only: no closure/discharge claims.
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

REPAIR_TEMPLATES = {
    0: "Gate-level fix: satisfy strict-side provenance/noncyclic/no-legacy/no-closure declarations.",
    1: "Probe-level fix (P1146): candidate sign/shape refinement on strict-side input family.",
    2: "Obstruction-level fix (P1147): address phase-induced first-flip mechanism with explicit strict premise.",
    3: "Phase-family fix (P1148): revise admissible phase-premise mapping and finite-domain constraints.",
    4: "Audit fix (P1149): restore reproducibility consistency across generated artifacts.",
}


def main() -> None:
    src = GEN / "p1153_strict_candidate_quality_ranking_summary.json"
    data = json.loads(src.read_text(encoding="utf-8"))

    plans = []
    for item in data.get("ranking", []):
        fs = item.get("failure_stage")
        idx = fs.get("index") if isinstance(fs, dict) else None
        if idx is None:
            action = "No repair required (currently admissible candidate-only)."
            priority = "monitor"
        else:
            action = REPAIR_TEMPLATES.get(idx, "Unknown stage: add explicit diagnostic mapping before repair.")
            priority = "high" if idx >= 2 else "medium"

        plans.append({
            "candidate": item.get("candidate"),
            "status": item.get("status"),
            "failure_stage_index": idx,
            "repair_action": action,
            "repair_priority": priority,
        })

    out = {
        "packet": "P1157",
        "as_of": "2026-05-10",
        "source_ranking": str(src),
        "repair_plans": plans,
        "note": "Methodological repair planning only; no strict-core closure or QW-2191 discharge claim.",
    }

    out_path = GEN / "p1157_failure_stage_repair_planner_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1157] planned {len(plans)} candidates, wrote {out_path}")


if __name__ == "__main__":
    main()
