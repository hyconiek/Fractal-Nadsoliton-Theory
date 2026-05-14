#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
IN1610 = GEN / "p1610_s560_full_strict_lagrangian_and_eom_chain_export_summary.json"
IN1613 = GEN / "p1613_s563_symbolic_euler_lagrange_identity_replay_summary.json"

def _load(p: Path) -> dict:
    if not p.exists():
        raise FileNotFoundError(p.name)
    return json.loads(p.read_text(encoding="utf-8"))

def main() -> None:
    s62 = _load(IN1562)
    s10 = _load(IN1610)
    s13 = _load(IN1613)

    closure = s13.get("strict_core_closure", {})
    ready = (
        s62.get("status", "").startswith("PASS")
        and s10.get("status", "").startswith("PASS")
        and s13.get("status", "").startswith("PASS")
        and not closure.get("missing_exports")
        and not closure.get("missing_witnesses")
        and not closure.get("missing_theorems")
    )

    dossier = {
        "kernel": s62.get("strict_kernel_params", {}),
        "coefficients": s62.get("derived_lagrangian_coefficients", {}),
        "full_lagrangian": s10.get("full_lagrangian_density", {}),
        "eom": s10.get("euler_lagrange_equations", {}),
        "identity_replay": s13.get("identity_replay", {}),
    }

    summary = {
        "checkpoint": "P1614_S564",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1614_STRICT_TOE_CLOSURE_DOSSIER" if ready else "KEEP_OPEN_P1614_DOSSIER_GAP",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "dossier": dossier,
        "strict_core_closure": {
            "status": "CLOSED" if ready else "OPEN",
            "missing_exports": closure.get("missing_exports", []),
            "missing_witnesses": closure.get("missing_witnesses", []),
            "missing_theorems": closure.get("missing_theorems", []),
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Attach reproducible symbolic notebook or CAS script as executable appendix to dossier.",
        "lay_summary": "To jeden kompletny dokument: od kernela strict aż po równania ruchu i sprawdzenie tożsamości, bez skrótów przez legacy.",
    }

    out = GEN / "p1614_s564_strict_toe_closure_dossier_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Wrote {out}")

if __name__ == "__main__":
    main()
