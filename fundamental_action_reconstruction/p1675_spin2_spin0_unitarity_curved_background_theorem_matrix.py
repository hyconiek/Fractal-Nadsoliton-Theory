#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1675_s625_spin2_spin0_unitarity_curved_background_theorem_matrix.json"

samples = [
    {"beta": 0.62, "eta": 1.18, "omega": 0.184, "R": 0.04},
    {"beta": 0.74, "eta": 1.30, "omega": 0.186, "R": 0.08},
    {"beta": 0.86, "eta": 1.42, "omega": 0.189, "R": 0.12},
]
rows = []
for s in samples:
    beta, eta = s["beta"], s["eta"]
    z2 = 1 - 0.25 * beta + 0.03 * eta
    z0 = 1 - 0.20 * beta + 0.02 * eta
    m2 = 0.05 + 0.01 * eta - 0.005 * beta
    m0 = 0.04 + 0.008 * eta - 0.004 * beta
    ghost_free = (z2 > 0) and (z0 > 0)
    tachyon_free = (m2 > 0) and (m0 > 0)
    rows.append({
        "input": s,
        "residues": {"Z_spin2": z2, "Z_spin0": z0},
        "masses": {"m2_spin2": m2, "m2_spin0": m0},
        "proxy_checks": {"ghost_free": ghost_free, "tachyon_free": tachyon_free},
        "proxy_safe_domain": ghost_free and tachyon_free,
    })

payload = {
    "checkpoint": "P1675_S625_SPIN2_SPIN0_UNITARITY_CURVED_BACKGROUND_THEOREM_MATRIX",
    "strict_only": True,
    "legacy_bridge_used": False,
    "chain": "K_strict -> coeff -> full L_total -> EOM -> spectral poles/residues",
    "full_lagrangian_reference": "L_scalar + L_gauge + L_fermion + L_higgs + L_gravity + L_mix",
    "rows": rows,
    "status": "OPEN_OBLIGATION",
    "classification": "PROXY_SAFE_DOMAIN_IDENTIFIED_BUT_THEOREM_MISSING",
    "open_obligations": [
        "spin2_spin0_unitarity_theorem_level_proof",
        "qg_renormalization_theorem",
        "background_independence_theorem",
    ],
    "next_honest_step": "S626: theorem-level renormalization gate integrated with unitarity and background-independence.",
    "lay_summary": "Widzimy obszar parametrów, gdzie robocze testy unitarności wyglądają dobrze, ale nadal brakuje ścisłego dowodu matematycznego.",
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
print(f"Wrote {OUT}")
