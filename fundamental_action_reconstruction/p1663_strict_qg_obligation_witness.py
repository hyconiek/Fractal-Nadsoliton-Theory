#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1663_s613_strict_qg_obligation_witness.json"

# strict kernel point (operational sample)
omega, phi, beta, eta, A = 0.18575, 4.5231, 0.72, 1.31, 1.0

coeff = {
    "Mpl2": A * (1.0 + beta),
    "cR2": A * beta / (1.0 + eta),
    "cRic2": A * beta * eta / (1.0 + eta),
    "cRiem2": A * beta * (1.0 + eta) / 4.0,
    "ZA": 1.0 + beta**2,
    "muH2": A * omega**2,
    "lambdaH": (1.0 + eta**2) / (1.0 + beta),
    "xiHR": beta / (1.0 + beta),
}

# 1-loop proxy counterterm map
counterterms = {
    "delta_cR2": 0.01 * coeff["cR2"],
    "delta_cRic2": 0.01 * coeff["cRic2"],
    "delta_cRiem2": 0.01 * coeff["cRiem2"],
}

# scalar profile and derivatives: h(x)=h0+h1*x+h2*x^2
h0, h1, h2 = 1.2, -0.4, 0.3
x = 0.7
h = h0 + h1 * x + h2 * x * x
hxx = 2.0 * h2

# EOM: h'' + muH2 h + lambdaH h^3 -2 xiHR R h = 0
# two backgrounds (proxy)
R_mink = 0.0
R_frw = 0.08

def scalar_residual(R: float) -> float:
    return hxx + coeff["muH2"] * h + coeff["lambdaH"] * h**3 - 2.0 * coeff["xiHR"] * R * h

res_mink = scalar_residual(R_mink)
res_frw = scalar_residual(R_frw)

# background-independence proxy: same operator structure, only R-value substitution
same_form = True

# spin-2 unitarity proxy (non-theorem)
unitarity_proxy_pass = (coeff["Mpl2"] > 0.0) and (coeff["cRic2"] >= 0.0)

payload = {
    "checkpoint": "P1663_S613_STRICT_QG_OBLIGATION_WITNESS",
    "strict_only": True,
    "legacy_bridge_used": False,
    "kernel_point": {"omega": omega, "phi": phi, "beta": beta, "eta": eta, "A": A},
    "coefficients": coeff,
    "counterterm_basis": ["R^2", "Ricci^2", "Riemann^2"],
    "counterterms_proxy": counterterms,
    "unitarity_proxy": {
        "condition": "Mpl2>0 and cRic2>=0",
        "pass": unitarity_proxy_pass,
    },
    "background_independence_proxy": {
        "backgrounds": ["Minkowski", "FRW_proxy"],
        "same_operator_form": same_form,
    },
    "eom_scalar_residuals": {
        "Minkowski_R0": res_mink,
        "FRW_R0p08": res_frw,
    },
    "status": "OPEN_OBLIGATION",
    "open_obligations": [
        "full_1loop_counterterm_derivation",
        "spin2_propagator_unitarity_theorem",
        "background_field_quantization_independence_proof",
    ],
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"Wrote {OUT}")
