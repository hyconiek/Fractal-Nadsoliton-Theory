#!/usr/bin/env python3
"""P1609/S559: strict full-chain closure audit from kernel to EOM and final strict-core closure."""
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"
IN1607 = GEN / "p1607_s557_final_g3_replay_with_np1_witness_bridge_summary.json"
IN1608 = GEN / "p1608_s558_frozen_strict_core_theorem_bundle_summary.json"


def _load(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    s62 = _load(IN1562)
    s63 = _load(IN1563)
    s07 = _load(IN1607)
    s08 = _load(IN1608)

    stages = {
        "kernel_to_coefficients": s62.get("status", "").startswith("PASS"),
        "coefficients_to_lagrangian_to_eom": s63.get("status", "").startswith("PASS"),
        "final_g3_replay": s07.get("status", "").startswith("PASS"),
        "frozen_bundle": s08.get("status", "").startswith("PASS"),
    }

    closure = s07.get("strict_core_closure", {})
    missing_exports = closure.get("missing_exports", [])
    missing_witnesses = closure.get("missing_witnesses", [])
    missing_theorems = closure.get("missing_theorems", [])

    closure_clean = not (missing_exports or missing_witnesses or missing_theorems)
    stage_ready = all(stages.values())
    strict_closed = stage_ready and closure_clean and not s08.get("bundle", {}).get("legacy_bridge_used", True)

    status = "PASS_P1609_STRICT_FULL_CHAIN_CLOSURE_AUDIT" if strict_closed else "KEEP_OPEN_P1609_STRICT_FULL_CHAIN_GAP"

    summary = {
        "checkpoint": "P1609_S559",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain_stages": stages,
        "strict_core_closure": {
            "status": "CLOSED" if strict_closed else "OPEN",
            "missing_exports": missing_exports,
            "missing_witnesses": missing_witnesses,
            "missing_theorems": missing_theorems,
        },
        "legacy_bridge_used": bool(s08.get("bundle", {}).get("legacy_bridge_used", False)),
        "external_team_validation_required": False,
        "next_honest_step": (
            "If CLOSED: continue with strict-only ToE refinement on this frozen baseline. "
            "If OPEN: discharge first missing element from exports/witnesses/theorems with theorem-level artifact."
        ),
        "lay_summary": (
            "To audyt pełnego łańcucha: od kernela strict do równań ruchu i bramki G3. "
            "Domknięcie jest uczciwe tylko jeśli nie ma żadnych brakujących eksportów, świadków ani twierdzeń."
        ),
    }

    out = GEN / "p1609_s559_strict_kernel_to_full_lagrangian_closure_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
