#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1664_s614_strict_full_lagrangian_manifest_and_inversion.json"

# kernel point
omega, phi, beta, eta, A = 0.18575, 4.5231, 0.72, 1.31, 1.0

coeff = {
    "Mpl2": A * (1.0 + beta),
    "cR2": A * beta / (1.0 + eta),
    "cRic2": A * beta * eta / (1.0 + eta),
    "cRiem2": A * beta * (1.0 + eta) / 4.0,
    "Z3": 1.0 + 0.8 * beta**2,
    "Z2": 1.0 + 0.6 * beta**2,
    "Z1": 1.0 + 0.4 * beta**2,
    "muH2": A * omega**2,
    "lambdaH": (1.0 + eta**2) / (1.0 + beta),
    "xiHR": beta / (1.0 + beta),
    "chiRG": beta * eta / (2.0 + beta),
    "yu": 0.9 * omega,
    "yd": 0.5 * omega,
    "ye": 0.3 * omega,
}

# inverse map (local)
beta_rec = coeff["xiHR"] / (1.0 - coeff["xiHR"])
A_rec = coeff["Mpl2"] / (1.0 + beta_rec)
omega_rec = math.sqrt(max(coeff["muH2"] / A_rec, 0.0))
eta_rec = math.sqrt(max(coeff["lambdaH"] * (1.0 + beta_rec) - 1.0, 0.0))

errs = {
    "beta_abs": abs(beta_rec - beta),
    "A_abs": abs(A_rec - A),
    "omega_abs": abs(omega_rec - omega),
    "eta_abs": abs(eta_rec - eta),
}

tol = 1e-12
local_pass = all(v < tol for v in errs.values())

lagrangian_manifest = {
    "L_GR": "sqrt(-g)[(Mpl2/2)R + cR2 R^2 + cRic2 Ricci^2 + cRiem2 Riemann^2]",
    "L_gauge": "-sqrt(-g)/4[Z3 G^2 + Z2 W^2 + Z1 B^2]",
    "L_H": "sqrt(-g)[(D H)^2 - muH2 H^†H - lambdaH(H^†H)^2]",
    "L_fermion": "sqrt(-g) sum_f i psibar_f gamma^a e_a^mu D_mu psi_f",
    "L_Yukawa": "-sqrt(-g)[yu Qbar H~ uR + yd Qbar H dR + ye Lbar H eR + h.c.]",
    "L_mix": "sqrt(-g)[xiHR(H^†H)R + chiRG R(G^2+W^2+B^2)]",
    "L_QG_ct": "sqrt(-g)[dcR2 R^2 + dcRic2 Ricci^2 + dcRiem2 Riemann^2]",
}

payload = {
    "checkpoint": "P1664_S614_STRICT_FULL_LAGRANGIAN_MANIFEST_AND_INVERSION",
    "strict_only": True,
    "legacy_bridge_used": False,
    "kernel_input": {"omega": omega, "phi": phi, "beta": beta, "eta": eta, "A": A},
    "coefficient_map": coeff,
    "lagrangian_manifest": lagrangian_manifest,
    "inverse_recovery": {
        "beta": beta_rec,
        "A": A_rec,
        "omega": omega_rec,
        "eta": eta_rec,
        "abs_errors": errs,
        "local_pass": local_pass,
        "tolerance": tol,
    },
    "status": "OPEN_OBLIGATION",
    "open_obligations": [
        "global_invertibility_theorem",
        "spin2_unitarity_operator_test",
        "full_1loop_counterterm_derivation",
        "background_independence_quantum_proof",
    ],
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"Wrote {OUT}")
print(f"local_pass={local_pass}")
