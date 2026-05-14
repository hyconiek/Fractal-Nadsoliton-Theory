#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1666_s616_curved_background_spin2_spin0_spectral_proxy.json"

omega, phi, beta, eta, A = 0.18575, 4.5231, 0.72, 1.31, 1.0
Mpl2 = A * (1 + beta)
cR2 = A * beta / (1 + eta)
cRic2 = A * beta * eta / (1 + eta)
cRiem2 = A * beta * (1 + eta) / 4

backgrounds = {
    "Minkowski": 0.0,
    "FRW_proxy": 0.08,
}

results = {}
for name, Rbg in backgrounds.items():
    den2 = cRic2 + 4 * cRiem2 + 0.1 * Rbg
    den0 = 6 * cR2 + 2 * cRic2 + 0.05 * Rbg
    m2 = Mpl2 / den2 if den2 != 0 else None
    m0 = Mpl2 / den0 if den0 != 0 else None
    results[name] = {
        "Rbg": Rbg,
        "m2_spin2": m2,
        "m2_spin0": m0,
        "tachyon_free_spin2": (m2 is not None and m2 > 0),
        "tachyon_free_spin0": (m0 is not None and m0 > 0),
        "residue_sign_spin2_proxy": -1 if den2 > 0 else 1,
        "residue_sign_spin0_proxy": 1 if den0 > 0 else -1,
    }

stable_signs = (
    results["Minkowski"]["residue_sign_spin2_proxy"] == results["FRW_proxy"]["residue_sign_spin2_proxy"]
    and results["Minkowski"]["residue_sign_spin0_proxy"] == results["FRW_proxy"]["residue_sign_spin0_proxy"]
)

payload = {
    "checkpoint": "P1666_S616_CURVED_BACKGROUND_SPIN2_SPIN0_SPECTRAL_PROXY",
    "strict_only": True,
    "legacy_bridge_used": False,
    "kernel_point": {"omega": omega, "phi": phi, "beta": beta, "eta": eta, "A": A},
    "coefficients": {"Mpl2": Mpl2, "cR2": cR2, "cRic2": cRic2, "cRiem2": cRiem2},
    "background_results": results,
    "stable_residue_signs_across_backgrounds": stable_signs,
    "status": "OPEN_OBLIGATION",
    "open_obligations": [
        "full_quadratic_kernel_diagonalization",
        "gauge_fixed_projector_residue_theorem",
        "global_unitarity_plus_renormalization_proof",
    ],
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
print(f"Wrote {OUT}")
