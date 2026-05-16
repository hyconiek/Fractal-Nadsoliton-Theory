#!/usr/bin/env python3
"""P1845 S795 strict full-Lagrangian explicit equation sheet checkpoint."""

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
    p1842 = load("p1842_s792_strict_kernel_to_full_lagrangian_formula_pack_checkpoint.json")
    p1839 = load("p1839_s789_strict_full_lagrangian_term_delivery_contract_checkpoint.json")

    pack = p1842.get("strict_formula_pack", {})
    sectors = p1839.get("term_delivery_contract", [])

    equation_sheet = {
        "kernel": pack.get("kernel_anchor", "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)"),
        "coefficient_map": pack.get("coefficient_map_schema", {}),
        "L_total_explicit_structure": (
            "L_total = L_GR[g,R,Ric,Riem] + L_gauge[F_A,F_W,F_G] + "
            "L_fermion[psi,D] + L_Higgs[H,DH,V(H)] + L_mix[H,chi] + "
            "L_int[gauge-fermion,gauge-higgs,yukawa,higher_mix] + "
            "L_covariant[sqrt(-g),tetrad,spin_connection]"
        ),
        "eom_generation_rule": "E_X := d/dx^mu(∂L_total/∂(∂_mu X)) - ∂L_total/∂X for each field block",
    }

    sector_table = [
        {
            "sector": s.get("sector", "UNKNOWN"),
            "density_symbol": s.get("density_symbol", "MISSING"),
            "missing_exports": s.get("missing_exports", []),
            "ready": s.get("sector_ready", False),
        }
        for s in sectors
    ]

    out = {
        "packet_id": "P1845",
        "stage_id": "S795",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "explicit_equation_sheet": equation_sheet,
        "sector_delivery_table": sector_table,
        "technical_progress": "Full-chain equation sheet is made explicit from K_strict through coefficient map and full L_total to EOM generation rule.",
        "proven": "The strict derivation grammar is explicit and sector-addressable, preventing skeleton-only ambiguity.",
        "open": "Sector term exports and theorem-level QG witnesses are still missing for closure.",
        "false_pass_risk": "Equation-sheet explicitness alone is not a theorem witness; missing sector exports keep closure open.",
        "next_honest_step": "Instantiate each sector equation term with explicit coefficients and export blockwise covariant EOM residual traces.",
        "lay_explanation": "To jawna karta równań: pokazuje jak przejść od kernela do pełnego wzoru i równań ruchu, ale brakuje jeszcze dowodów końcowych.",
    }

    path = GEN / "p1845_s795_strict_full_lagrangian_explicit_equation_sheet_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
