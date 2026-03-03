#!/usr/bin/env python3
"""
QW-1716: Unified uncertainty budget for headline FIN observables.

Domains:
1) Electroweak angle
2) Fine-structure inverse
3) Mass sector (selected particles)

Uncertainty decomposition:
sigma_tot^2 = sigma_num^2 + sigma_model^2 + sigma_calib^2 + sigma_data^2
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1716_uncertainty_budget.json"
OUT_MD = ROOT / "RAPORT_QW1716_UNCERTAINTY_BUDGET.md"


def sqrt_sum(*vals: float) -> float:
    return math.sqrt(sum(v * v for v in vals))


def rel_pct(abs_err: float, ref: float) -> float:
    if ref == 0:
        return float("nan")
    return abs(abs_err) / abs(ref) * 100.0


def main() -> None:
    # Base parameters
    alpha_geo = 4.0 * math.log(2.0)
    beta_tors = 0.01
    gamma = 1.52
    m_top = 173_000.0

    # Experimental references / assumed data uncertainties
    sin2_exp = 0.23122
    sin2_data_sigma = 0.00003

    alpha_inv_exp = 137.035999
    alpha_inv_data_sigma = 0.000021

    masses_exp = {
        "Top": (173_000.0, 300.0),      # MeV, conservative sigma
        "Bottom": (4_180.0, 30.0),
        "Tau": (1_776.9, 0.2),
        "Charm": (1_270.0, 20.0),
        "Muon": (105.7, 0.0002),
        "Electron": (0.511, 1e-6),
    }
    q_map = {"Top": 0, "Bottom": 7, "Tau": 9, "Charm": 9, "Muon": 14, "Electron": 24}

    # Effective corrections from QW-1708
    delta_w = +7.398420630e-04
    delta_vac = -1.509312219e-03

    # --- Electroweak observable ---
    sin2_pred = (alpha_geo / 12.0) * (1.0 + delta_w)
    sin2_res = sin2_pred - sin2_exp
    # component uncertainties
    sigma_num_sin2 = 1e-10
    sigma_model_sin2 = abs(sin2_res) * 0.15   # residual model uncertainty share
    sigma_calib_sin2 = abs((alpha_geo / 12.0) * delta_w) * 0.20
    sigma_data_sin2 = sin2_data_sigma
    sigma_tot_sin2 = sqrt_sum(sigma_num_sin2, sigma_model_sin2, sigma_calib_sin2, sigma_data_sin2)

    # --- Fine structure observable ---
    alpha_inv_pred = (alpha_geo / (2.0 * beta_tors)) * (1.0 - beta_tors) * (1.0 + delta_vac)
    alpha_inv_res = alpha_inv_pred - alpha_inv_exp
    sigma_num_alpha = 1e-9
    sigma_model_alpha = abs(alpha_inv_res) * 0.15
    sigma_calib_alpha = abs((alpha_geo / (2.0 * beta_tors)) * (1.0 - beta_tors) * delta_vac) * 0.20
    sigma_data_alpha = alpha_inv_data_sigma
    sigma_tot_alpha = sqrt_sum(sigma_num_alpha, sigma_model_alpha, sigma_calib_alpha, sigma_data_alpha)

    # --- Mass sector ---
    # Use Delta model from QW-1711 as effective correction profile:
    l1, l2, l3 = -0.096670, -0.031035, +0.023793
    sector = {"Top": 0.0, "Bottom": 1.0, "Tau": -1.0, "Charm": 1.0, "Muon": -1.0, "Electron": -1.0}
    gen = {"Top": 3.0, "Bottom": 3.0, "Tau": 3.0, "Charm": 2.0, "Muon": 2.0, "Electron": 1.0}

    mass_rows = []
    for p, (exp_val, data_sigma) in masses_exp.items():
        q = q_map[p]
        base = m_top * (4.0 ** (-(gamma * q / 4.0)))
        delta = l1 * (q / 24.0) + l2 * sector[p] + l3 * (gen[p] - 2.0)
        pred = base * math.exp(delta)
        res = pred - exp_val

        sigma_num = max(1e-10 * exp_val, 1e-8)
        sigma_model = abs(res) * 0.50  # dominant unknown mechanism share
        sigma_calib = abs(pred - base) * 0.25
        sigma_data = data_sigma
        sigma_tot = sqrt_sum(sigma_num, sigma_model, sigma_calib, sigma_data)

        mass_rows.append(
            {
                "particle": p,
                "pred_mev": pred,
                "exp_mev": exp_val,
                "residual_mev": res,
                "residual_pct": rel_pct(res, exp_val),
                "sigma_num_mev": sigma_num,
                "sigma_model_mev": sigma_model,
                "sigma_calib_mev": sigma_calib,
                "sigma_data_mev": sigma_data,
                "sigma_tot_mev": sigma_tot,
            }
        )

    # Dominance map
    def dominant_component(row: dict) -> str:
        comps = {
            "num": row["sigma_num_mev"],
            "model": row["sigma_model_mev"],
            "calib": row["sigma_calib_mev"],
            "data": row["sigma_data_mev"],
        }
        return max(comps, key=comps.get)

    for r in mass_rows:
        r["dominant_component"] = dominant_component(r)

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_params": {
            "alpha_geo": alpha_geo,
            "beta_tors": beta_tors,
            "gamma": gamma,
            "m_top_mev": m_top,
        },
        "effective_corrections": {"delta_w": delta_w, "delta_vac": delta_vac, "l1": l1, "l2": l2, "l3": l3},
        "electroweak": {
            "sin2_pred": sin2_pred,
            "sin2_exp": sin2_exp,
            "residual": sin2_res,
            "sigma_num": sigma_num_sin2,
            "sigma_model": sigma_model_sin2,
            "sigma_calib": sigma_calib_sin2,
            "sigma_data": sigma_data_sin2,
            "sigma_tot": sigma_tot_sin2,
        },
        "alpha_inverse": {
            "pred": alpha_inv_pred,
            "exp": alpha_inv_exp,
            "residual": alpha_inv_res,
            "sigma_num": sigma_num_alpha,
            "sigma_model": sigma_model_alpha,
            "sigma_calib": sigma_calib_alpha,
            "sigma_data": sigma_data_alpha,
            "sigma_tot": sigma_tot_alpha,
        },
        "mass_sector": mass_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1716: UNCERTAINTY BUDGET",
        "",
        f"- Data UTC: {output['generated_utc']}",
        "",
        "## 1) Electroweak",
        f"- sin²θ residual = {sin2_res:+.6e}",
        f"- sigma_tot = {sigma_tot_sin2:.6e}",
        "",
        "## 2) Fine structure inverse",
        f"- residual = {alpha_inv_res:+.6e}",
        f"- sigma_tot = {sigma_tot_alpha:.6e}",
        "",
        "## 3) Mass sector (residual %, dominant sigma)",
    ]
    for r in mass_rows:
        lines.append(
            f"- {r['particle']}: residual={r['residual_pct']:+.2f}%, sigma_tot={r['sigma_tot_mev']:.4f} MeV, dominant={r['dominant_component']}"
        )
    lines.extend(
        [
            "",
            "## Wniosek",
            "- W sektorze mas dominuje niepewność modelowa/kalibracyjna, nie numeryczna.",
            "- Dalsze domknięcie wymaga poprawy mechanizmu flavor/mass, nie tylko lepszej numerki.",
            "",
            "## Artefakty",
            f"- JSON szczegółowy: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1716] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1716] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
