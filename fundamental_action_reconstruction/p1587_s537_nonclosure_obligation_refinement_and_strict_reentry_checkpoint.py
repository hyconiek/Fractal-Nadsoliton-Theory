#!/usr/bin/env python3
"""P1587/S537: refine nonclosure obligations and strict re-entry plan from strict kernel chain."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1586 = GEN / "p1586_s536_final_toe_composition_theorem_or_nonclosure_certificate_summary.json"


def main() -> None:
    if not IN1586.exists():
        raise FileNotFoundError("Run P1586 before P1587")

    s = json.loads(IN1586.read_text(encoding="utf-8"))
    final_obj = s["final_object"]
    is_open = final_obj["status"] == "OPEN"

    if is_open:
        missing = final_obj["nonclosure_certificate"]["remaining_missing"]
        prioritized = [
            "G1_full_domain_selector_gap_proof",
            "G2_global_SM_GR_stability_proof",
            "G3_final_ToE_composition_export",
        ]
        status = "KEEP_OPEN_P1587_STRICT_REENTRY_PLAN_EXPORTED"
    else:
        missing = []
        prioritized = []
        status = "PASS_P1587_NO_REENTRY_NEEDED"

    summary = {
        "checkpoint": "P1587_S537",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": s["strict_chain"],
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "nonclosure_refinement": {
            "open": is_open,
            "inherited_missing": missing,
            "prioritized_obligations": prioritized,
            "reentry_policy": "No legacy bridge; return to strict proof objects only",
        },
        "strict_core_closure": {
            "status": "OPEN" if is_open else "CLOSED",
            "qw2191_closed": False if is_open else True,
            "toe_closed": False if is_open else True,
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1588_discharge_G1_selector_gap_on_full_domain_then_compose_G2",
        "lay_summary": "Jeśli ToE nie jest domknięte, checkpoint porządkuje brakujące dowody i ustawia uczciwy plan powrotu bez skrótów i bez legacy bridge."
    }

    out = GEN / "p1587_s537_nonclosure_obligation_refinement_and_strict_reentry_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
