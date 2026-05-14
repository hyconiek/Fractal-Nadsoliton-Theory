#!/usr/bin/env python3
"""P1586/S536: final strict ToE composition theorem export or nonclosure certificate."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1585 = GEN / "p1585_s535_formal_l1_l3_theorem_objects_and_toe_composition_summary.json"


def main() -> None:
    if not IN1585.exists():
        raise FileNotFoundError("Run P1585 before P1586")

    s = json.loads(IN1585.read_text(encoding="utf-8"))
    gate_ready = bool(s["toe_composition_gate"]["ready"])

    if gate_ready:
        status = "PASS_P1586_FINAL_STRICT_TOE_COMPOSITION_THEOREM_EXPORTED"
        closure = {
            "status": "CLOSED",
            "qw2191_closed": True,
            "toe_closed": True,
            "export": "T1586_final_strict_ToE_composition_theorem"
        }
    else:
        status = "KEEP_OPEN_P1586_NONCLOSURE_CERTIFICATE"
        closure = {
            "status": "OPEN",
            "qw2191_closed": False,
            "toe_closed": False,
            "nonclosure_certificate": {
                "reason": "Formal L1/L3 global obligations not fully discharged",
                "remaining_missing": s["strict_core_closure"]["remaining_missing"],
            }
        }

    summary = {
        "checkpoint": "P1586_S536",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": s["strict_chain"],
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "final_object": closure,
        "external_team_validation_required": False,
        "next_honest_step": "If OPEN: construct missing global proofs; if CLOSED: independent internal replay.",
        "lay_summary": "Wydano finalny theorem tylko jeśli bramka formalna była gotowa; w przeciwnym razie wydano certyfikat uczciwego niedomknięcia."
    }

    out = GEN / "p1586_s536_final_toe_composition_theorem_or_nonclosure_certificate_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
