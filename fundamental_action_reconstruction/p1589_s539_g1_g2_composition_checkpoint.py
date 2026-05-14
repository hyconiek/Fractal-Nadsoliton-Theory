#!/usr/bin/env python3
"""P1589/S539: compose G1 selector-gap discharge with G2 global stability candidate."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1588 = GEN / "p1588_s538_g1_full_domain_selector_gap_discharge_summary.json"
IN1584 = GEN / "p1584_s534_l1_witness_refinement_and_l3_global_stability_object_summary.json"


def main() -> None:
    if not IN1588.exists() or not IN1584.exists():
        raise FileNotFoundError("Run P1588 and P1584 before P1589")

    g1 = json.loads(IN1588.read_text(encoding="utf-8"))
    g2 = json.loads(IN1584.read_text(encoding="utf-8"))

    g1_pass = bool(g1["G1_selector_gap"]["pass"])
    g2_pass = bool(g2["l3_global_stability_object_candidate"]["pass"])
    composed = g1_pass and g2_pass

    status = "PASS_P1589_G1_G2_COMPOSITION_READY" if composed else "KEEP_OPEN_P1589_G1_G2_NOT_READY"

    summary = {
        "checkpoint": "P1589_S539",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": g1["strict_chain"],
        "composition": {
            "G1_selector_gap_pass": g1_pass,
            "G2_global_stability_pass": g2_pass,
            "ready_for_final_export": composed,
        },
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_core_closure": {
            "status": "OPEN",
            "qw2191_closed": False,
            "remaining_missing": [] if composed else ["G1/G2 full theorem-level strengthening", "G3 final ToE export"],
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1590_export_G3_final_ToE_object_if_ready_else_issue_refined_nonclosure",
        "lay_summary": "Połączono dwa kluczowe warunki (G1 i G2), aby ocenić gotowość do finalnego eksportu ToE bez skrótów metodologicznych."
    }

    out = GEN / "p1589_s539_g1_g2_composition_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
