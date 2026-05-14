#!/usr/bin/env python3
"""P1607/S557: final strict G3 replay using W1606/T1606 imports and closure/nonclosure decision."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1606 = GEN / "p1606_s556_np1_witness_and_bridge_upgrade_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"


def main() -> None:
    if not IN1606.exists() or not IN1563.exists():
        raise FileNotFoundError("Run P1606 and P1563 before P1607")

    s06 = json.loads(IN1606.read_text(encoding="utf-8"))
    s63 = json.loads(IN1563.read_text(encoding="utf-8"))

    witness_ready = bool(s06.get("np1_witness_object", {}).get("exported", False))
    bridge_ready = bool(s06.get("bridge_upgrade_object", {}).get("exported", False))
    eom_ready = s63.get("status", "").startswith("PASS")

    g3_ready = witness_ready and bridge_ready and eom_ready
    status = "PASS_P1607_FINAL_STRICT_CLOSURE_READY" if g3_ready else "KEEP_OPEN_P1607_FINAL_STRICT_NONCLOSURE"

    summary = {
        "checkpoint": "P1607_S557",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s06.get("strict_chain", {}),
        "final_replay_gate": {
            "w1606_import_ready": witness_ready,
            "t1606_import_ready": bridge_ready,
            "eom_ready": eom_ready,
            "g3_ready": g3_ready,
        },
        "strict_core_closure": {
            "status": "CLOSED" if g3_ready else "OPEN",
            "missing_exports": [] if eom_ready else ["E_lagrangian_to_eom_export"],
            "missing_witnesses": [] if witness_ready else ["W_G1_full_domain_selector_gap_discharge"],
            "missing_theorems": [] if bridge_ready and g3_ready else ["T_G3_final_strict_ToE_composition"],
        },
        "external_team_validation_required": False,
        "next_honest_step": "If CLOSED: emit frozen strict-core theorem bundle; if OPEN: reopen only failing gate with explicit nonbridge theorem obligations.",
        "lay_summary": "To finalna replay-bramka: jeśli świadek NP1, upgrade bridge i równania ruchu są gotowe, możemy domknąć strict-core; inaczej zostaje uczciwe nonclosure."
    }

    out = GEN / "p1607_s557_final_g3_replay_with_np1_witness_bridge_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
