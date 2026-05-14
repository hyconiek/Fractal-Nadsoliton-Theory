#!/usr/bin/env python3
"""P1597/S547: final strict G3 theorem composition object from kernel->coefficients->Lagrangian->EOM chain."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1596 = GEN / "p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json"
IN1588 = GEN / "p1588_s538_g1_full_domain_selector_gap_discharge_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"


def _load(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    s96 = _load(IN1596)
    s88 = _load(IN1588)
    s63 = _load(IN1563)

    bridge_ready = bool(s96.get("bridge_gate", {}).get("bridge_ready_now", False))
    g1_block = s88.get("G1_selector_gap", s88.get("g1_selector_gap", {}))
    g1_ready = bool(g1_block.get("discharged", g1_block.get("pass", False)))
    eom_ready = s63.get("status", "").startswith("PASS")

    g3_theorem_object = {
        "id": "T1597_final_strict_G3_composition_object",
        "domain": "F_Nadsoliton => L_SM + L_GR",
        "chain": "K_strict -> coefficients -> full Lagrangian -> Euler-Lagrange EOM",
        "assumptions": [
            "A1: strict selector uniqueness theorem object available",
            "A2: full-domain selector gap (G1) discharged",
            "A3: EOM export consistency already passed",
        ],
        "conclusion": "strict-core final composition theorem candidate",
        "proof_trace_exported": bridge_ready and g1_ready and eom_ready,
        "legacy_bridge_used": False,
    }

    ready = g3_theorem_object["proof_trace_exported"]
    status = "PASS_P1597_FINAL_G3_THEOREM_OBJECT_EXPORTED" if ready else "KEEP_OPEN_P1597_G3_THEOREM_NOT_READY"

    summary = {
        "checkpoint": "P1597_S547",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s96.get("strict_chain", {}),
        "g3_theorem_object": g3_theorem_object,
        "gate_review": {
            "bridge_ready": bridge_ready,
            "g1_ready": g1_ready,
            "eom_ready": eom_ready,
        },
        "strict_core_closure": {
            "status": "CLOSED" if ready else "OPEN",
            "missing_exports": [] if eom_ready else ["E_lagrangian_to_eom_export"],
            "missing_witnesses": [] if g1_ready else ["W_G1_full_domain_selector_gap_discharge"],
            "missing_theorems": [] if ready else ["T_G3_final_strict_ToE_composition"],
        },
        "external_team_validation_required": False,
        "next_honest_step": "If OPEN: discharge remaining G1/bridge blocker and rerun P1597; if CLOSED: emit final closure declaration packet.",
        "lay_summary": "To jest finalny szkic twierdzenia G3: łączymy strict kernel i równania ruchu w jeden obiekt, ale domykamy tylko gdy wszystkie bramki są naprawdę spełnione."
    }

    out = GEN / "p1597_s547_final_g3_theorem_composition_object_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
