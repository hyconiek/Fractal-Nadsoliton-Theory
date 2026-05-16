#!/usr/bin/env python3
"""P1839 S789 strict full-Lagrangian term-delivery contract checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1836 = load("p1836_s786_strict_full_lagrangian_non_skeleton_manifest_checkpoint.json")
    p1834 = load("p1834_s784_strict_full_lagrangian_sector_export_gate_checkpoint.json")
    p1838 = load("p1838_s788_strict_kernel_to_eom_to_qg_theorem_gate_map_checkpoint.json")

    manifest = p1836.get("full_lagrangian_non_skeleton_manifest", {}).get("L_total", {})
    sector_gate = {s.get("sector", "UNKNOWN"): s for s in p1834.get("full_lagrangian_sector_gate", [])}

    term_delivery = []
    for sector, payload in manifest.items():
        gate = sector_gate.get(sector, {})
        required = gate.get("required_exports", [])
        term_delivery.append(
            {
                "sector": sector,
                "density_symbol": payload.get("density_symbol", "MISSING_DENSITY_SYMBOL"),
                "required_exports": required,
                "delivered_exports": gate.get("delivered_exports", []),
                "missing_exports": gate.get("missing_exports", required),
                "sector_ready": gate.get("sector_gate_pass", False),
            }
        )

    tg_map = p1838.get("theorem_gate_map", {})

    out = {
        "packet_id": "P1839",
        "stage_id": "S789",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "term_delivery_contract": term_delivery,
        "theorem_gate_dependencies": tg_map,
        "technical_progress": "Full-Lagrangian sector manifest is now transformed into a term-delivery contract linked to TG dependency map.",
        "proven": "TG gates cannot be satisfied unless sector-level term-delivery contract is fulfilled.",
        "open": "Term-level deliveries remain incomplete across sectors, so theorem gates stay open.",
        "false_pass_risk": "Declaring TG readiness without completing sector term-delivery contract would bypass strict physics content requirements.",
        "next_honest_step": "Complete missing term exports in gravity/gauge/fermion first, then propagate updates to TG dependency map re-evaluation.",
        "lay_explanation": "Każdy dział wzoru ma listę brakujących składników; dopiero po jej zamknięciu można uczciwie przejść do bramek teorii.",
    }

    path = GEN / "p1839_s789_strict_full_lagrangian_term_delivery_contract_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
