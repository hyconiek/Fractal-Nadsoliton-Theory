#!/usr/bin/env python3
"""P1595/S545: strict-only final G3 composition gate on kernel→coefficients→Lagrangian→EOM route."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1593 = GEN / "p1593_s543_focused_first_obligation_discharge_summary.json"
IN1594 = GEN / "p1594_s544_focused_second_obligation_g2_discharge_summary.json"
IN1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"
IN1582 = GEN / "p1582_s532_strict_selector_uniqueness_theorem_bridge_to_full_lagrangian_summary.json"


def _load(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    s93 = _load(IN1593)
    s94 = _load(IN1594)
    s62 = _load(IN1562)
    s63 = _load(IN1563)
    s82 = _load(IN1582)

    g1_ready = bool(s93["focused_obligation"]["pass"])
    g2_ready = bool(s94["focused_obligation"]["pass"])
    coeff_export_ready = s62.get("status", "").startswith("PASS")
    eom_export_ready = s63.get("status", "").startswith("PASS")
    theorem_bridge_ready = bool(s82.get("strict_uniqueness_bridge", {}).get("bridge_ready", False))

    ready = all([g1_ready, g2_ready, coeff_export_ready, eom_export_ready, theorem_bridge_ready])

    missing_exports = []
    if not coeff_export_ready:
        missing_exports.append("E_kernel_to_coefficients_export")
    if not eom_export_ready:
        missing_exports.append("E_lagrangian_to_eom_export")

    missing_witnesses = []
    if not g1_ready:
        missing_witnesses.append("W_G1_full_domain_selector_gap_discharge")
    if not g2_ready:
        missing_witnesses.append("W_G2_global_stability_object")

    missing_theorems = []
    if not theorem_bridge_ready:
        missing_theorems.append("T_selector_uniqueness_to_full_lagrangian_bridge")
    if not ready:
        missing_theorems.append("T_G3_final_strict_ToE_composition")

    status = "PASS_P1595_G3_FINAL_EXPORT_CANDIDATE" if ready else "KEEP_OPEN_P1595_G3_NOT_READY"

    summary = {
        "checkpoint": "P1595_S545",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": {
            "kernel": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
            "coefficients": s62.get("coefficients", {}),
            "lagrangian": s62.get("lagrangian", {}),
            "eom": s63.get("eom_export", {}),
        },
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only, no legacy bridge)",
        "preconditions": {
            "G1_ready": g1_ready,
            "G2_ready": g2_ready,
            "coeff_export_ready": coeff_export_ready,
            "eom_export_ready": eom_export_ready,
            "theorem_bridge_ready": theorem_bridge_ready,
        },
        "G3": {
            "candidate_export_ready": ready,
            "object": "T1595_final_strict_ToE_composition_candidate" if ready else None,
        },
        "strict_core_closure": {
            "status": "CLOSED" if ready else "OPEN",
            "missing_exports": missing_exports,
            "missing_witnesses": missing_witnesses,
            "missing_theorems": missing_theorems,
        },
        "external_team_validation_required": False,
        "next_honest_step": "Discharge first missing item from missing_theorems with explicit theorem object and proof trace.",
        "lay_summary": "To jest końcowa bramka strict: pokazujemy cały tor od kernela do równań ruchu i uczciwie wypisujemy, czego brakuje do pełnego domknięcia ToE."
    }

    out = GEN / "p1595_s545_final_g3_attempt_from_g1_g2_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
