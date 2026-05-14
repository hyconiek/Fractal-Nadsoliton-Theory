#!/usr/bin/env python3
"""P1591/S541: internal replay package if closed, or focused gap-work package if open."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1590 = GEN / "p1590_s540_g3_final_export_or_refined_nonclosure_summary.json"


def main() -> None:
    if not IN1590.exists():
        raise FileNotFoundError("Run P1590 before P1591")

    s = json.loads(IN1590.read_text(encoding="utf-8"))
    closed = s["strict_core_closure"]["status"] == "CLOSED"

    if closed:
        status = "PASS_P1591_INTERNAL_REPLAY_PACKAGE_READY"
        package = {
            "mode": "internal_replay",
            "artifacts": [
                "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json",
                "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json",
                "p1590_s540_g3_final_export_or_refined_nonclosure_summary.json",
            ],
        }
    else:
        status = "KEEP_OPEN_P1591_FOCUSED_GAP_WORK_PACKAGE"
        package = {
            "mode": "gap_focus",
            "targets": s["strict_core_closure"].get("remaining_missing", []),
            "priority": ["G1", "G2", "G3"],
        }

    summary = {
        "checkpoint": "P1591_S541",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": s["strict_chain"],
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "work_package": package,
        "strict_core_closure": s["strict_core_closure"],
        "external_team_validation_required": False,
        "next_honest_step": "If OPEN: execute focused theorem discharge; if CLOSED: run deterministic replay and freeze package.",
        "lay_summary": "Checkpoint przekłada status ToE na konkretny pakiet pracy: albo domykanie braków, albo replikację gotowego wyniku."
    }

    out = GEN / "p1591_s541_internal_replay_or_gap_focus_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
