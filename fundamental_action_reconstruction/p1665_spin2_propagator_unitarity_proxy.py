#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1665_s615_spin2_propagator_unitarity_proxy.json"

# strict sample from prior chain
omega, phi, beta, eta, A = 0.18575, 4.5231, 0.72, 1.31, 1.0
Mpl2 = A * (1.0 + beta)
cR2 = A * beta / (1.0 + eta)
cRic2 = A * beta * eta / (1.0 + eta)
cRiem2 = A * beta * (1.0 + eta) / 4.0

# proxy spin-2 effective denominator coefficient
den_spin2 = cRic2 + 4.0 * cRiem2
m2_spin2 = Mpl2 / den_spin2 if den_spin2 != 0.0 else None

massless_ok = Mpl2 > 0.0
massive_nontachyonic = (m2_spin2 is not None) and (m2_spin2 > 0.0)

# sign proxy for heavy residue (true => possible ghost risk)
ghost_risk_flag = den_spin2 > 0.0

payload = {
    "checkpoint": "P1665_S615_SPIN2_PROPAGATOR_UNITARITY_PROXY",
    "strict_only": True,
    "legacy_bridge_used": False,
    "kernel_point": {"omega": omega, "phi": phi, "beta": beta, "eta": eta, "A": A},
    "gr_coefficients": {
        "Mpl2": Mpl2,
        "cR2": cR2,
        "cRic2": cRic2,
        "cRiem2": cRiem2,
    },
    "spin2_proxy": {
        "denominator_coeff": den_spin2,
        "m2_spin2": m2_spin2,
        "massless_positive_energy": massless_ok,
        "massive_nontachyonic": massive_nontachyonic,
        "ghost_risk_flag": ghost_risk_flag,
    },
    "status": "OPEN_OBLIGATION",
    "open_obligations": [
        "full_operator_propagator_derivation",
        "residue_sign_theorem_on_curved_backgrounds",
        "complete_unitarity_proof_spin2_spin0",
    ],
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"Wrote {OUT}")
