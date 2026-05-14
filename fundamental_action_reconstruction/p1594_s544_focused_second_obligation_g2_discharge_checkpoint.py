#!/usr/bin/env python3
"""P1594/S544: focused discharge attempt for second obligation (G2 global stability object)."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1593 = GEN / "p1593_s543_focused_first_obligation_discharge_summary.json"
IN1584 = GEN / "p1584_s534_l1_witness_refinement_and_l3_global_stability_object_summary.json"


def main() -> None:
    if not IN1593.exists() or not IN1584.exists():
        raise FileNotFoundError("Run P1593 and P1584 before P1594")

    s93 = json.loads(IN1593.read_text(encoding="utf-8"))
    s84 = json.loads(IN1584.read_text(encoding="utf-8"))

    c = s84["l3_global_stability_object_candidate"]
    g2_pass = bool(c["pass"]) and float(c["coercivity_margin"]) > 0.11

    remaining = s93["strict_core_closure"].get("remaining_missing", [])
    target = remaining[0] if remaining else "none"

    status = "PASS_P1594_G2_DISCHARGED_CANDIDATE" if g2_pass else "KEEP_OPEN_P1594_G2_NOT_DISCHARGED"

    summary = {
        "checkpoint": "P1594_S544",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": s93["strict_chain"],
        "focused_obligation": {
            "target": target,
            "discharge_candidate": "G2_global_SM_GR_stability_theorem_object_candidate",
            "metrics": c,
            "pass": g2_pass,
        },
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_core_closure": {
            "status": "OPEN",
            "remaining_missing": [m for m in remaining if m != target] if g2_pass else remaining,
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1595: attempt final G3 composition export using discharged G1/G2 candidates.",
        "lay_summary": "Checkpoint atakuje drugi główny brak: globalny obiekt stabilności. Jeśli przejdzie, zostaje już tylko finalna kompozycja ToE."
    }

    out = GEN / "p1594_s544_focused_second_obligation_g2_discharge_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
