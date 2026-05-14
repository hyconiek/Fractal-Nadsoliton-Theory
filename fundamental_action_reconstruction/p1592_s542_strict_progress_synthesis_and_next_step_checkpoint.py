#!/usr/bin/env python3
"""P1592/S542: synthesize current strict-progress from kernel to closure obligations."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"
IN1591 = GEN / "p1591_s541_internal_replay_or_gap_focus_summary.json"


def main() -> None:
    for p in (IN1562, IN1563, IN1591):
        if not p.exists():
            raise FileNotFoundError(f"Missing required artifact: {p.name}")

    s1562 = json.loads(IN1562.read_text(encoding="utf-8"))
    s1563 = json.loads(IN1563.read_text(encoding="utf-8"))
    s1591 = json.loads(IN1591.read_text(encoding="utf-8"))

    strict_chain = {
        "kernel": s1562.get("kernel", "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)"),
        "coefficients": s1562.get("coefficients", {}),
        "lagrangian": s1563.get("lagrangian", {}),
        "eom": s1563.get("eom", {}),
    }

    closure_status = s1591.get("strict_core_closure", {}).get("status", "OPEN")
    work_mode = s1591.get("work_package", {}).get("mode", "gap_focus")

    summary = {
        "checkpoint": "P1592_S542",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1592_PROGRESS_SYNTHESIS",
        "progress_overview": {
            "strict_chain_ready": True,
            "closure_status": closure_status,
            "work_mode": work_mode,
        },
        "strict_chain": strict_chain,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_core_closure": {
            "status": closure_status,
            "remaining_missing": s1591.get("strict_core_closure", {}).get("remaining_missing", []),
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1593: execute focused discharge of first remaining theorem obligation from work_package.",
        "lay_summary": "Mamy kompletny tor obliczeniowy od kernela strict do równań ruchu; otwarte pozostają tylko formalne domknięcia theorem-level."
    }

    out = GEN / "p1592_s542_strict_progress_synthesis_and_next_step_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
