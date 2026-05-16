#!/usr/bin/env python3
"""P1822 S772 strict S1 evidence-pack readiness checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def is_witness_ready(status: str) -> bool:
    return status.startswith("PASS_ZERO")


def main() -> None:
    p1820 = load("p1820_s770_strict_current_priority_bottleneck_to_execution_contract_checkpoint.json")
    p1821 = load("p1821_s771_strict_priority_lane_sequencer_checkpoint.json")

    snapshot = p1820.get("current_priority_snapshot", {})
    lane = p1821.get("execution_lane", [])
    s1 = next((x for x in lane if x.get("step") == "S1_EXPORT_EA_EH_ELG_BOUNDARY_H1_BW_INPUTS"), {})

    s1_targets = s1.get("targets", [])

    target_status = {t: snapshot.get(t, "OPEN_UNKNOWN") for t in s1_targets}
    ready_flags = {t: is_witness_ready(st) for t, st in target_status.items()}

    missing = [t for t, ok in ready_flags.items() if not ok]

    out = {
        "packet_id": "P1822",
        "stage_id": "S772",
        "status": "PASS_ZERO" if not missing else "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "lane_step": "S1_EXPORT_EA_EH_ELG_BOUNDARY_H1_BW_INPUTS",
        "targets": s1_targets,
        "target_status": target_status,
        "witness_readiness": ready_flags,
        "missing_witness_targets": missing,
        "completion_rule": s1.get("completion_rule", "all targets must have residual/witness traces attached"),
        "technical_progress": "S1 step converted to machine-checkable readiness gate for TG1 input packaging.",
        "proven": "Readiness can be evaluated per target without theorem-gate promotion.",
        "open": "At least one S1 target still lacks PASS_ZERO-level witness trace attachment.",
        "false_pass_risk": "Skipping per-target readiness allows TG1 run with partial evidence and ambiguous failures.",
        "next_honest_step": "Attach residual/witness traces for each missing S1 target, then rerun P1822 before TG1 execution.",
        "lay_explanation": "To checklista przed głównym testem: każdy wymagany element musi mieć własny dowód, inaczej nie idziemy dalej.",
    }

    path = GEN / "p1822_s772_strict_s1_evidence_pack_readiness_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
