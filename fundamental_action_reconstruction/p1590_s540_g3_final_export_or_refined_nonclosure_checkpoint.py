#!/usr/bin/env python3
"""P1590/S540: export G3 final ToE object if G1&G2 ready, else refined nonclosure certificate."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1589 = GEN / "p1589_s539_g1_g2_composition_summary.json"


def main() -> None:
    if not IN1589.exists():
        raise FileNotFoundError("Run P1589 before P1590")

    s = json.loads(IN1589.read_text(encoding="utf-8"))
    ready = bool(s["composition"]["ready_for_final_export"])

    if ready:
        status = "PASS_P1590_G3_FINAL_TOE_OBJECT_EXPORTED"
        g3 = {
            "type": "T1590_final_strict_ToE_export",
            "toe_closed": True,
            "qw2191_closed": True,
        }
        remaining = []
    else:
        status = "KEEP_OPEN_P1590_REFINED_NONCLOSURE"
        g3 = {
            "type": "T1590_nonclosure_certificate_refined",
            "toe_closed": False,
            "qw2191_closed": False,
            "reason": "G1/G2 not jointly theorem-ready",
        }
        remaining = s["strict_core_closure"]["remaining_missing"]

    summary = {
        "checkpoint": "P1590_S540",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": s["strict_chain"],
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "G3_object": g3,
        "strict_core_closure": {
            "status": "CLOSED" if ready else "OPEN",
            "remaining_missing": remaining,
        },
        "external_team_validation_required": False,
        "next_honest_step": "If OPEN: refine G1/G2 theorem objects; if CLOSED: internal replay package.",
        "lay_summary": "Finalny obiekt ToE jest eksportowany tylko przy pełnej gotowości G1 i G2; w przeciwnym razie system uczciwie raportuje niedomknięcie."
    }

    out = GEN / "p1590_s540_g3_final_export_or_refined_nonclosure_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
