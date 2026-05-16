#!/usr/bin/env python3
"""P1834 S784 strict full-Lagrangian sector export gate checkpoint."""

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
    p1833 = load("p1833_s783_strict_full_lagrangian_and_reverse_closure_worklist_checkpoint.json")

    sectors = p1833.get("full_lagrangian_sector_worklist", [])
    sector_gate = []
    for s in sectors:
        required = s.get("required_exports", [])
        sector_gate.append(
            {
                "sector": s.get("sector", "UNKNOWN"),
                "required_exports": required,
                "delivered_exports": [],
                "missing_exports": required,
                "sector_gate_pass": False,
            }
        )

    pass_all = all(x["sector_gate_pass"] for x in sector_gate) and len(sector_gate) > 0

    out = {
        "packet_id": "P1834",
        "stage_id": "S784",
        "status": "PASS_ZERO" if pass_all else "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "full_lagrangian_sector_gate": sector_gate,
        "gate_rule": "all sectors must pass before promoting forward chain to theorem-gate input",
        "reverse_chain_unlock_rule": "forward full-lagrangian sector gate PASS is required before reverse Helmholtz/BRST/Cutkosky closure attempts",
        "technical_progress": "Sector-level full-Lagrangian export obligations are now enforced by an explicit gate object.",
        "proven": "Forward full-Lagrangian chain cannot be treated as complete until each sector passes export gate checks.",
        "open": "All sectors currently fail due to missing explicit export attachments.",
        "false_pass_risk": "Using aggregate forward labels without per-sector export gate can hide missing physics terms.",
        "next_honest_step": "Attach explicit exports sector-by-sector (starting gravity/gauge/fermion) and rerun P1834 before reverse theorem gates.",
        "lay_explanation": "Każdy sektor pełnego Lagrangianu musi zdać własną kontrolę kompletności, zanim wolno iść dalej.",
    }

    path = GEN / "p1834_s784_strict_full_lagrangian_sector_export_gate_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
