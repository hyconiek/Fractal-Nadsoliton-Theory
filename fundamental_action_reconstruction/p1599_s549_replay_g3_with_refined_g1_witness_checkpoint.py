#!/usr/bin/env python3
"""P1599/S549: replay final G3 gate using refined strict G1 witness and strict bridge status."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1596 = GEN / "p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json"
IN1598 = GEN / "p1598_s548_strict_g1_witness_refinement_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"


def _load(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    s96 = _load(IN1596)
    s98 = _load(IN1598)
    s63 = _load(IN1563)

    bridge_ready = bool(s96.get("bridge_gate", {}).get("bridge_ready_now", False))
    g1_ready = bool(s98.get("g1_refined_witness", {}).get("pass", False))
    eom_ready = s63.get("status", "").startswith("PASS")

    g3_ready = bridge_ready and g1_ready and eom_ready
    status = "PASS_P1599_REPLAY_G3_READY_FOR_CLOSURE" if g3_ready else "KEEP_OPEN_P1599_REPLAY_G3_NOT_READY"

    summary = {
        "checkpoint": "P1599_S549",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s96.get("strict_chain", {}),
        "replay_gate": {
            "bridge_ready": bridge_ready,
            "g1_refined_witness_ready": g1_ready,
            "eom_ready": eom_ready,
            "g3_ready": g3_ready,
        },
        "strict_core_closure": {
            "status": "CLOSED" if g3_ready else "OPEN",
            "missing_exports": [] if eom_ready else ["E_lagrangian_to_eom_export"],
            "missing_witnesses": [] if g1_ready else ["W_G1_full_domain_selector_gap_discharge"],
            "missing_theorems": [] if g3_ready else ["T_G3_final_strict_ToE_composition"],
        },
        "external_team_validation_required": False,
        "next_honest_step": "If OPEN: target first failing gate in replay_gate and export its dedicated strict artifact; if CLOSED: emit final strict-core closure declaration.",
        "lay_summary": "To jest powtórka końcowej bramki G3 po wzmocnieniu świadka G1: jeśli wszystkie bramki są zielone, można domykać; jeśli nie, wiemy dokładnie co blokuje ToE."
    }

    out = GEN / "p1599_s549_replay_g3_with_refined_g1_witness_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
