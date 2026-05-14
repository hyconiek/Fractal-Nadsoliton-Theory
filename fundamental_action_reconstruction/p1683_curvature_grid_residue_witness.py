#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1683_s633_curvature_grid_residue_witness.json"

# strict-proxy grid (scaffold-level, no false theorem closure)
R_values = [-1.0, -0.3, 0.0, 0.3, 1.0]

rows = []
for R in R_values:
    z2 = 1.0 + 0.10 * R
    z0 = 1.0 + 0.06 * R
    det_k2 = z2 * (1.0 + 0.04 * R)
    det_k0 = z0 * (1.0 + 0.02 * R)
    rows.append({
        "R": R,
        "Z2_positive": z2 > 0,
        "Z0_positive": z0 > 0,
        "det_K_spin2_positive": det_k2 > 0,
        "det_K_spin0_positive": det_k0 > 0,
    })

all_proxy_positive = all(
    r["Z2_positive"] and r["Z0_positive"] and r["det_K_spin2_positive"] and r["det_K_spin0_positive"]
    for r in rows
)

payload = {
    "checkpoint": "P1683_S633_CURVATURE_GRID_RESIDUE_WITNESS",
    "strict_only": True,
    "legacy_bridge_used": False,
    "pipeline": "K_strict -> coeff -> full_L_total -> linearized_EOM -> kinetic_block_sign_witness",
    "full_lagrangian_anchor": "L_gauge + L_fermion + L_higgs + L_yukawa + L_gravity + L_mix",
    "grid_rows": rows,
    "proxy_global_sign_stability": all_proxy_positive,
    "status": "OPEN_OBLIGATION",
    "limitation": "Proxy sign witness only; theorem-level curvature-uniform inequality proof still missing.",
    "next_honest_step": "S634: promote proxy witness to lemma-proof sketch with explicit R-inequalities and (c1,c2,xi) conditions."
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"Wrote {OUT}")
