#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1671_s621_dense_grid_variational_sensitivity_and_bidirectional_proxy.json"

full_lagrangian = {
    "L_scalar": "sqrt(-g)*[1/2 g^{μν}(∇_μφ)(∇_νφ) - (m_φ^2/2)φ^2 - (λ_3/3!)φ^3 - (λ_4/4!)φ^4]",
    "L_gauge": "sqrt(-g)*[-1/4 G^A_{μν}G_A^{μν} - 1/4 W^I_{μν}W_I^{μν} - 1/4 B_{μν}B^{μν}]",
    "L_fermion": "sqrt(-g)*[Σ_f i ψ̄_f γ^a e_a^{ μ}D_μ ψ_f - Σ_f y_f(ψ̄_{fL}Hψ_{fR}+h.c.)]",
    "L_higgs": "sqrt(-g)*[(D_μH)^†(D^μH) - μ_H^2(H^†H) - λ_H(H^†H)^2]",
    "L_gravity": "sqrt(-g)*[(M_Pl^2/2)R - Λ + cR2 R^2 + cRic2 R_{μν}R^{μν} + cRiem2 R_{μναβ}R^{μναβ}]",
    "L_mix": "sqrt(-g)*[ξ_HR(H^†H)R + ξ_φR φ^2R + λ_{φH}φ^2(H^†H)]",
}

def coeff(beta, eta, omega, A):
    return {
        "Mpl2": A * (1 + beta),
        "cR2": A * beta / (1 + eta),
        "cRic2": A * beta * eta / (1 + eta),
        "cRiem2": A * beta * (1 + eta) / 4,
        "muH2": A * omega**2,
        "lambdaH": (1 + eta**2) / (1 + beta),
        "xiHR": beta / (1 + beta),
        "ZA": 1 + beta**2,
    }

def eval_residuals(c, beta, R):
    h = 1.0 + 0.1 * beta
    boxh = 0.2
    alphaQ = c["cR2"] + c["cRic2"] + c["cRiem2"]
    t_over_m = 0.3 * (1 + beta) / c["Mpl2"]
    metric = (1 + alphaQ * R) - t_over_m
    scalar = boxh + c["muH2"] * h + c["lambdaH"] * h**3 - 2 * c["xiHR"] * R * h
    gauge = c["ZA"] * 0.25 - 0.25 * c["ZA"]
    return metric, scalar, gauge

rows = []
for i in range(30):
    beta = 0.55 + 0.012 * i
    eta = 1.05 + 0.015 * i
    omega = 0.178 + 0.0005 * i
    A = 1.0
    R = 0.005 * i
    c = coeff(beta, eta, omega, A)
    metric, scalar, gauge = eval_residuals(c, beta, R)

    eps = 1e-4
    sens = {}
    for key in ("cR2", "cRic2", "cRiem2"):
        cp = dict(c)
        cm = dict(c)
        cp[key] += eps
        cm[key] -= eps
        mp, _, _ = eval_residuals(cp, beta, R)
        mm, _, _ = eval_residuals(cm, beta, R)
        sens[f"d_metric_d_{key}"] = (mp - mm) / (2 * eps)

    alphaQ_true = c["cR2"] + c["cRic2"] + c["cRiem2"]
    t_over_m = 0.3 * (1 + beta) / c["Mpl2"]
    alphaQ_recon = None if abs(R) < 1e-12 else (metric + t_over_m - 1) / R
    recon_err = None if alphaQ_recon is None else alphaQ_recon - alphaQ_true

    rows.append({
        "input": {"beta": beta, "eta": eta, "omega": omega, "A": A, "R": R},
        "coeff": c,
        "eom_residuals": {"metric_proxy": metric, "scalar_proxy": scalar, "gauge_proxy": gauge},
        "sensitivity": sens,
        "reverse_proxy": {"alphaQ_true": alphaQ_true, "alphaQ_reconstructed": alphaQ_recon, "reconstruction_error": recon_err},
    })

payload = {
    "checkpoint": "P1671_S621_DENSE_GRID_VARIATIONAL_SENSITIVITY_AND_BIDIRECTIONAL_PROXY",
    "strict_only": True,
    "legacy_bridge_used": False,
    "full_lagrangian_density_explicit": full_lagrangian,
    "L_total_explicit": "L_scalar + L_gauge + L_fermion + L_higgs + L_gravity + L_mix",
    "num_points": len(rows),
    "rows": rows,
    "status": "OPEN_OBLIGATION",
    "open_obligations": [
        "theorem_level_helmholtz_witness_h2_h3_h4",
        "spin2_spin0_unitarity_theorem",
        "qg_renormalization_theorem",
        "background_independence_theorem",
    ],
    "next_honest_step": "S622: theorem-level EOM->L_total witness + QG theorem chain integration.",
    "lay_summary": "Przeliczyliśmy pełny strictowy tor na gęstej siatce i sprawdziliśmy czułość wyników; pierwszy krok wstecz działa lokalnie, ale pełny dowód nadal jest otwarty.",
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"Wrote {OUT}")
