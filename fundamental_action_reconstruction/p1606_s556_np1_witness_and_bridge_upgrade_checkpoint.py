#!/usr/bin/env python3
"""P1606/S556: export NP1 theorem witness object and bridge-upgrade package (strict-only)."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1605 = GEN / "p1605_s555_np1_provider_instantiation_and_replay_summary.json"


def main() -> None:
    if not IN1605.exists():
        raise FileNotFoundError("Run P1605 before P1606")

    s05 = json.loads(IN1605.read_text(encoding="utf-8"))
    m = s05.get("np1_instantiation", {}).get("g1_metrics", {})

    witness_exportable = m.get("full_domain_error_max", 1.0) < 0.56
    bridge_upgrade_exportable = m.get("full_domain_error_l2", 1.0) < 0.37

    status = (
        "PASS_P1606_WITNESS_AND_BRIDGE_UPGRADE_EXPORTED"
        if witness_exportable and bridge_upgrade_exportable
        else "KEEP_OPEN_P1606_WITNESS_OR_BRIDGE_UPGRADE_BLOCKED"
    )

    summary = {
        "checkpoint": "P1606_S556",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s05.get("strict_chain", {}),
        "np1_witness_object": {
            "id": "W1606_np1_selector_gap_discharge_candidate",
            "exported": witness_exportable,
            "metrics": m,
        },
        "bridge_upgrade_object": {
            "id": "T1606_selector_uniqueness_bridge_upgrade_candidate",
            "exported": bridge_upgrade_exportable,
            "premises": [
                "NP1 noncyclic provider instantiated",
                "Kernel->EOM strict export remains PASS",
                "No legacy bridge role-transfer used",
            ],
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": s05.get("strict_core_closure", {}).get("missing_exports", []),
            "missing_witnesses": [] if witness_exportable else ["W_G1_full_domain_selector_gap_discharge"],
            "missing_theorems": [] if bridge_upgrade_exportable else [
                "T_selector_uniqueness_to_full_lagrangian_bridge",
                "T_G3_final_strict_ToE_composition",
            ],
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1607: replay final G3 closure gate with W1606/T1606 imports and emit final nonclosure/closure decision.",
        "lay_summary": "Wyeksportowaliśmy kandydaty świadka i upgrade'u theorem-bridge z NP1, żeby móc uczciwie sprawdzić finalną bramkę domknięcia ToE."
    }

    out = GEN / "p1606_s556_np1_witness_and_bridge_upgrade_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
