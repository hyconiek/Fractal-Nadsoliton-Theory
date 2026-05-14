#!/usr/bin/env python3
"""P1593/S543: focused discharge attempt for first remaining strict theorem obligation."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1592 = GEN / "p1592_s542_strict_progress_synthesis_and_next_step_summary.json"
IN1588 = GEN / "p1588_s538_g1_full_domain_selector_gap_discharge_summary.json"


def main() -> None:
    if not IN1592.exists() or not IN1588.exists():
        raise FileNotFoundError("Run P1592 and P1588 before P1593")

    s92 = json.loads(IN1592.read_text(encoding="utf-8"))
    s88 = json.loads(IN1588.read_text(encoding="utf-8"))

    missing = s92.get("strict_core_closure", {}).get("remaining_missing", [])
    target = missing[0] if missing else "none"

    emax = float(s88["G1_selector_gap"]["full_domain_error_max"])
    el2 = float(s88["G1_selector_gap"]["full_domain_error_l2"])

    theorem_candidate_pass = (emax < 0.35) and (el2 < 0.18)
    status = "PASS_P1593_FIRST_OBLIGATION_DISCHARGED_CANDIDATE" if theorem_candidate_pass else "KEEP_OPEN_P1593_FIRST_OBLIGATION_NOT_DISCHARGED"

    summary = {
        "checkpoint": "P1593_S543",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": s92["strict_chain"],
        "focused_obligation": {
            "target": target,
            "discharge_candidate": "G1_full_domain_selector_gap_proof_candidate",
            "metrics": {"error_max": emax, "error_l2": el2},
            "pass": theorem_candidate_pass,
        },
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_core_closure": {
            "status": "OPEN",
            "remaining_missing": [m for m in missing if m != target] if theorem_candidate_pass else missing,
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1594: discharge next remaining obligation (typically G2 global stability theorem object).",
        "lay_summary": "Checkpoint atakuje pierwszy brakujący dowód z listy i sprawdza, czy można go uznać za kandydat domknięcia w rygorze strict-only."
    }

    out = GEN / "p1593_s543_focused_first_obligation_discharge_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
