#!/usr/bin/env python3
"""
QW-2025: Mass-operator reformulation scan on eta-kernel branch.

Goal:
- after QW-2024 (flavor rework), test whether mass branch can be repaired
  with one shared extended mass operator (global coefficients only).
"""

from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2025_mass_operator_reformulation_scan_on_eta_kernel.json"
OUT_MD = ROOT / "RAPORT_QW2025_MASS_OPERATOR_REFORMULATION_SCAN_ON_ETA_KERNEL.md"


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def mass_eval(kernel: Dict[str, float], gamma_scale: float, coeffs: Dict[str, float]) -> Dict[str, object]:
    particles = [
        ("Top", 0.0, 173_000.0, 0.0, 3.0),
        ("Bottom", 7.0, 4_180.0, 1.0, 3.0),
        ("Tau", 9.0, 1_776.9, -1.0, 3.0),
        ("Charm", 9.0, 1_270.0, 1.0, 2.0),
        ("Muon", 14.0, 105.7, -1.0, 2.0),
        ("Electron", 24.0, 0.511, -1.0, 1.0),
    ]

    d1, d4 = np.array([1.0]), np.array([4.0])
    k1 = float(abs(kernel_fn(d1, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])[0]))
    k4 = float(abs(kernel_fn(d4, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])[0]))
    gamma_kernel = float(-4.0 * math.log(max(k4 / max(k1, 1e-15), 1e-15), 4.0) / 3.0)
    gamma_eff = float(np.clip(gamma_scale * gamma_kernel, 0.40, 2.60))

    rows: List[Dict[str, float]] = []
    errs = []
    for name, q, m_exp, sector, gen in particles:
        xq = q / 24.0
        xg = gen - 2.0

        base = 173_000.0 * (4.0 ** (-(gamma_eff * q / 4.0)))
        delta = (
            coeffs["c_q"] * xq
            + coeffs["c_s"] * sector
            + coeffs["c_g"] * xg
            + coeffs["c_q2"] * (xq**2)
            + coeffs["c_sg"] * sector * xg
        )

        pred = float(base * math.exp(delta))
        err = float(abs(pred - m_exp) / max(m_exp, 1e-15) * 100.0)
        errs.append(err)
        rows.append({
            "particle": name,
            "exp_mev": m_exp,
            "pred_mev": pred,
            "rel_err_pct": err,
        })

    return {
        "rows": rows,
        "gamma_kernel": gamma_kernel,
        "gamma_eff": gamma_eff,
        "mean_rel_err_pct": float(np.mean(errs)),
        "max_rel_err_pct": float(np.max(errs)),
    }


def main() -> None:
    d2021 = json.loads((ROOT / "report_qw2021_v2_eta_operator_beta_constraint_scan.json").read_text(encoding="utf-8"))
    d2024 = json.loads((ROOT / "report_qw2024_eta_kernel_isospin_flavor_rework_scan.json").read_text(encoding="utf-8"))

    kernel = d2021["selected"]["fit"]
    best2024 = d2024["best_row"]

    thresholds = {
        "mass_mean_rel_pct_max": 15.0,
        "mass_max_rel_pct_max": 75.0,
    }

    gamma_grid = np.linspace(0.45, 1.05, 7)
    coeff_vals = [-0.8, -0.4, 0.0, 0.4, 0.8]

    best = None
    pass_count = 0
    n_eval = 0

    for gamma_scale in gamma_grid:
        for c_q, c_s, c_g, c_q2, c_sg in itertools.product(coeff_vals, repeat=5):
            coeffs = {
                "c_q": float(c_q),
                "c_s": float(c_s),
                "c_g": float(c_g),
                "c_q2": float(c_q2),
                "c_sg": float(c_sg),
            }
            m = mass_eval(kernel=kernel, gamma_scale=float(gamma_scale), coeffs=coeffs)
            n_eval += 1

            flags = {
                "mass_mean_rel_pct_le_max": bool(m["mean_rel_err_pct"] <= thresholds["mass_mean_rel_pct_max"]),
                "mass_max_rel_pct_le_max": bool(m["max_rel_err_pct"] <= thresholds["mass_max_rel_pct_max"]),
            }
            all_pass = bool(all(flags.values()))
            if all_pass:
                pass_count += 1

            score = (
                m["mean_rel_err_pct"] / thresholds["mass_mean_rel_pct_max"]
                + 0.35 * (m["max_rel_err_pct"] / thresholds["mass_max_rel_pct_max"])
                + 0.05 * abs(coeffs["c_q2"])  # mild complexity penalty
                + 0.05 * abs(coeffs["c_sg"])
            )

            row = {
                "gamma_scale": float(gamma_scale),
                "coeffs": coeffs,
                "mass": m,
                "flags": flags,
                "all_pass": all_pass,
                "score": float(score),
            }
            if best is None or row["score"] < best["score"]:
                best = row

    verdict = "MASS_OPERATOR_REFORMULATION_PASS" if pass_count > 0 else "MASS_OPERATOR_REFORMULATION_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2021_v2_eta_operator_beta_constraint_scan.json:selected.fit",
        "kernel": kernel,
        "context_from_qw2024_best_row": {
            "p_amp": best2024["p_amp"],
            "r_dist": best2024["r_dist"],
            "flavor_ckm_mean_rel_pct": best2024["flavor"]["ckm_mean_rel_pct"],
            "flavor_pmns_mean_rel_pct": best2024["flavor"]["pmns_mean_rel_pct"],
            "gw_control_median_gap": best2024["gw"]["control_median_gap"],
        },
        "thresholds": thresholds,
        "search": {
            "n_eval": int(n_eval),
            "gamma_grid_size": int(len(gamma_grid)),
            "coeff_values": coeff_vals,
        },
        "pass_count_all_flags": int(pass_count),
        "best_row": best,
        "verdict": verdict,
        "required_next_step": (
            "COMBINE_WITH_FLAVOR_REWORK_IN_JOINT_SCAN"
            if pass_count > 0
            else "ESCALATE_TO_JOINT_MASS_FLAVOR_OPERATOR_REFORMULATION"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    b = best
    lines = [
        "# RAPORT QW-2025: MASS OPERATOR REFORMULATION SCAN ON ETA KERNEL",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- n_eval: {n_eval}",
        f"- pass_count_all_flags: {pass_count}",
        "",
        "## Best Row",
        f"- gamma_scale: {b['gamma_scale']:.4f}",
        f"- coeffs: {b['coeffs']}",
        f"- mass mean/max rel%: {b['mass']['mean_rel_err_pct']:.3f}/{b['mass']['max_rel_err_pct']:.3f}",
        f"- all_pass: {b['all_pass']}",
        "",
        "## QW-2024 Context (fixed branch)",
        f"- p_amp/r_dist: {best2024['p_amp']:.3f}/{best2024['r_dist']:.3f}",
        f"- flavor CKM/PMNS mean rel%: {best2024['flavor']['ckm_mean_rel_pct']:.3f}/{best2024['flavor']['pmns_mean_rel_pct']:.3f}",
        f"- GW control gap: {best2024['gw']['control_median_gap']:.6f}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2025] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2025] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2025] verdict={verdict} pass_count={pass_count}/{n_eval}")


if __name__ == "__main__":
    main()
