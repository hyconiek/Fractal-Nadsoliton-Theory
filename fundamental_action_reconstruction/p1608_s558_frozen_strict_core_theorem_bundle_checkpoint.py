#!/usr/bin/env python3
"""P1608/S558: freeze strict-core theorem bundle after final G3 replay decision."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"
IN1606 = GEN / "p1606_s556_np1_witness_and_bridge_upgrade_summary.json"
IN1607 = GEN / "p1607_s557_final_g3_replay_with_np1_witness_bridge_summary.json"


def _load(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    s62 = _load(IN1562)
    s63 = _load(IN1563)
    s06 = _load(IN1606)
    s07 = _load(IN1607)

    closed = s07.get("strict_core_closure", {}).get("status") == "CLOSED"

    bundle = {
        "bundle_id": "B1608_frozen_strict_core_theorem_bundle",
        "frozen_at_utc": datetime.now(UTC).isoformat(),
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "artifacts": {
            "kernel_to_coefficients": IN1562.name,
            "coefficients_to_lagrangian_to_eom": IN1563.name,
            "np1_witness_bridge_upgrade": IN1606.name,
            "final_g3_replay": IN1607.name,
        },
        "closure_status": "CLOSED" if closed else "OPEN",
        "legacy_bridge_used": False,
    }

    status = "PASS_P1608_STRICT_CORE_BUNDLE_FROZEN" if closed else "KEEP_OPEN_P1608_BUNDLE_NOT_FREEZABLE"

    summary = {
        "checkpoint": "P1608_S558",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "bundle": bundle,
        "strict_core_closure": s07.get("strict_core_closure", {}),
        "external_team_validation_required": False,
        "next_honest_step": "If CLOSED: maintain frozen baseline and open only post-closure extensions; if OPEN: return to failing gate with explicit theorem obligation.",
        "lay_summary": "To formalne zamrożenie wersji dowodowej strict-core: zapisujemy komplet artefaktów, które doprowadziły do aktualnego werdyktu domknięcia."
    }

    out = GEN / "p1608_s558_frozen_strict_core_theorem_bundle_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
