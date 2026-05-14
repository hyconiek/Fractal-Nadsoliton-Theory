#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1667_s617_gauge_fixed_quadratic_kernel_diagonalization_proxy.json"

omega, phi, beta, eta, A = 0.18575, 4.5231, 0.72, 1.31, 1.0
Mpl2 = A * (1 + beta)
cR2 = A * beta / (1 + eta)
cRic2 = A * beta * eta / (1 + eta)
cRiem2 = A * beta * (1 + eta) / 4

# proxy quadratic kernel block K = [[a,b],[b,d]]
a = cRic2 + 4*cRiem2
d = 6*cR2 + 2*cRic2
b = 0.15 * (cR2 + cRic2)

tr = a + d
det = a*d - b*b
rad = max(tr*tr - 4*det, 0.0)
l1 = 0.5 * (tr + math.sqrt(rad))
l2 = 0.5 * (tr - math.sqrt(rad))

m2_1 = Mpl2 / l1 if l1 != 0 else None
m2_2 = Mpl2 / l2 if l2 != 0 else None

payload = {
  "checkpoint": "P1667_S617_GAUGE_FIXED_QUADRATIC_KERNEL_DIAGONALIZATION_PROXY",
  "strict_only": True,
  "legacy_bridge_used": False,
  "kernel_point": {"omega": omega, "phi": phi, "beta": beta, "eta": eta, "A": A},
  "quadratic_kernel_proxy": {"a": a, "b": b, "d": d, "trace": tr, "det": det},
  "eigenvalues": {"lambda1": l1, "lambda2": l2},
  "mass_scales_proxy": {"m2_1": m2_1, "m2_2": m2_2},
  "nontachyonic_proxy": {"mode1": (m2_1 is not None and m2_1 > 0), "mode2": (m2_2 is not None and m2_2 > 0)},
  "residue_sign_proxy": {"mode1": -1 if l1 > 0 else 1, "mode2": -1 if l2 > 0 else 1},
  "status": "OPEN_OBLIGATION",
  "open_obligations": [
    "full_projector_level_diagonalization",
    "curved_background_residue_theorem",
    "combined_unitarity_renormalization_proof"
  ]
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
print(f"Wrote {OUT}")
