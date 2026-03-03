#!/usr/bin/env python3
"""
QW-1756: Omega-beta decoupling residual test.

Method:
- Start from QW-1752 good runs.
- Remove beta-envelope effect using beta posterior from QW-1755.
- Re-estimate omega from residualized spatial phase signal.
- Check whether omega sensitivity to beta perturbation drops.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1756_omega_beta_decoupling_residual_test.json"
OUT_MD = ROOT / "RAPORT_QW1756_OMEGA_BETA_DECOUPLING_RESIDUAL_TEST.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def analytic_signal_real(y: np.ndarray) -> np.ndarray:
    n = len(y)
    Y = np.fft.fft(y)
    H = np.zeros(n)
    if n % 2 == 0:
        H[0] = 1
        H[n // 2] = 1
        H[1 : n // 2] = 2
    else:
        H[0] = 1
        H[1 : (n + 1) // 2] = 2
    z = np.fft.ifft(Y * H)
    return z


def fit_phase_slope_hilbert(y_real: np.ndarray) -> Tuple[float, float, float]:
    d = np.arange(1, len(y_real) + 1, dtype=float)
    z = analytic_signal_real(y_real)
    amp = np.abs(z)
    ph = np.unwrap(np.angle(z))
    m = amp >= np.quantile(amp, 0.35)
    if np.sum(m) < 4:
        m = amp >= np.quantile(amp, 0.2)
    x = d[m]
    y = ph[m]
    A = np.column_stack([x, np.ones_like(x)])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope, intercept = float(coef[0]), float(coef[1])
    rmse = float(np.sqrt(np.mean((y - (slope * x + intercept)) ** 2)))
    return slope, intercept, rmse


def rankdata(x: np.ndarray) -> np.ndarray:
    idx = np.argsort(x)
    r = np.empty_like(idx, dtype=float)
    r[idx] = np.arange(len(x), dtype=float)
    return r


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    rx = rankdata(x)
    ry = rankdata(y)
    c = np.corrcoef(rx, ry)
    return float(c[0, 1])


def main() -> None:
    d1752 = load("report_qw1752_omega_orthogonal_observable.json")
    d1755 = load("report_qw1755_beta_null_vs_positive_evidence.json")

    beta_med = float(d1755["weighted_beta_posterior"]["median"])
    beta_ci = d1755["weighted_beta_posterior"]["ci95"]
    b_lo = float(beta_ci[0])
    b_hi = float(beta_ci[1])

    p = d1752.get("protocol", {}).get("good_filter", {})
    smin = float(p.get("snr_like_min", 4.0))
    rcmax = float(p.get("phase_rmse_complex_max", 0.45))
    rhmax = float(p.get("phase_rmse_hilbert_max", 0.55))

    rows = d1752.get("rows", [])
    good = [
        r
        for r in rows
        if np.isfinite(float(r.get("omega_hat", float("nan"))))
        and float(r.get("snr_like", 0.0)) >= smin
        and float(r.get("phase_rmse_complex", 1e9)) <= rcmax
        and float(r.get("phase_rmse_hilbert", 1e9)) <= rhmax
    ]
    if len(good) < 8:
        raise RuntimeError("Too few good rows from QW-1752 for residual decoupling test.")

    beta_grid = np.array(sorted(set([b_lo, beta_med, b_hi])), dtype=float)
    beta_grid = np.clip(beta_grid, 1e-6, 0.3)

    out_rows: List[Dict[str, object]] = []
    for r in good:
        yr = r.get("response_real", {})
        y = np.array([float(yr[str(i)]) for i in range(1, 13)], dtype=float)
        if np.max(np.abs(y)) > 1e-12:
            y = y / np.max(np.abs(y))
        d = np.arange(1, len(y) + 1, dtype=float)

        omega_grid = []
        rmse_grid = []
        for b in beta_grid:
            ycorr = y * (1.0 + b * d)
            s, ph, rm = fit_phase_slope_hilbert(ycorr)
            omega_grid.append(abs(s))
            rmse_grid.append(rm)

        omega_grid = np.array(omega_grid, dtype=float)
        rmse_grid = np.array(rmse_grid, dtype=float)
        i_ref = int(np.argmin(np.abs(beta_grid - beta_med)))
        omega_ref = float(omega_grid[i_ref])
        sens = float(np.max(np.abs(omega_grid - omega_ref)))

        out_rows.append(
            {
                "seed": int(r.get("seed", -1)),
                "n": int(r.get("n", -1)),
                "omega_hat_baseline": float(r.get("omega_hat")),
                "omega_hat_residual_beta_med": omega_ref,
                "omega_shift_abs": float(abs(omega_ref - float(r.get("omega_hat")))),
                "omega_beta_sensitivity_residual": sens,
                "phase_rmse_residual_median": float(np.median(rmse_grid)),
                "omega_grid": {f"{b:.8f}": float(o) for b, o in zip(beta_grid, omega_grid)},
            }
        )

    o0 = np.array([float(r["omega_hat_baseline"]) for r in out_rows], dtype=float)
    oc = np.array([float(r["omega_hat_residual_beta_med"]) for r in out_rows], dtype=float)
    os = np.array([float(r["omega_beta_sensitivity_residual"]) for r in out_rows], dtype=float)
    osh = np.array([float(r["omega_shift_abs"]) for r in out_rows], dtype=float)
    prm = np.array([float(r["phase_rmse_residual_median"]) for r in out_rows], dtype=float)

    omega_ci = [float(np.quantile(oc, 0.025)), float(np.quantile(oc, 0.975))]
    omega_iqr = float(np.quantile(oc, 0.75) - np.quantile(oc, 0.25))
    sens_q90 = float(np.quantile(os, 0.9))
    shift_med = float(np.median(osh))
    rmse_med = float(np.median(prm))
    rho = spearman(oc, os)

    pass_n = len(out_rows) >= 10
    pass_sens = sens_q90 <= 0.06
    pass_corr = np.isfinite(rho) and abs(rho) <= 0.35
    pass_spread = omega_iqr <= 0.30
    pass_nonboundary = omega_ci[0] > 0.10 and omega_ci[1] < 1.6
    pass_fit = rmse_med <= 0.40

    if pass_n and pass_sens and pass_corr and pass_spread and pass_nonboundary and pass_fit:
        verdict = "OMEGA_BETA_DECOUPLING_SUPPORTED"
    elif pass_n and pass_spread and pass_fit and (pass_sens or pass_corr):
        verdict = "OMEGA_BETA_DECOUPLING_PARTIAL"
    else:
        verdict = "OMEGA_BETA_DECOUPLING_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "n_rows_1752_good": len(out_rows),
            "beta_median_from_1755": beta_med,
            "beta_ci95_from_1755": [b_lo, b_hi],
            "beta_grid_used": [float(b) for b in beta_grid],
        },
        "summary": {
            "omega_residual_ci95": omega_ci,
            "omega_residual_iqr": omega_iqr,
            "omega_beta_sensitivity_q90": sens_q90,
            "omega_shift_abs_median_vs_1752": shift_med,
            "median_phase_rmse_residual": rmse_med,
            "spearman_omega_vs_sensitivity": rho,
            "baseline_omega_median_1752": float(np.median(o0)),
            "residual_omega_median": float(np.median(oc)),
        },
        "pass_flags": {
            "n_good_enough": bool(pass_n),
            "sensitivity_control": bool(pass_sens),
            "orthogonality_proxy_corr": bool(pass_corr),
            "spread_control": bool(pass_spread),
            "nonboundary_omega": bool(pass_nonboundary),
            "fit_quality": bool(pass_fit),
        },
        "verdict": verdict,
        "rows": out_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1756: OMEGA BETA DECOUPLING RESIDUAL TEST",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Good rows: {len(out_rows)}",
        f"- beta med/CI95 from 1755: {beta_med:.6f} / [{b_lo:.6f}, {b_hi:.6f}]",
        f"- residual omega CI95: [{omega_ci[0]:.6f}, {omega_ci[1]:.6f}]",
        f"- residual omega IQR: {omega_iqr:.6f}",
        f"- sensitivity q90: {sens_q90:.6f}",
        f"- spearman(omega,sensitivity): {rho:.6f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- n_good_enough: {pass_n}",
        f"- sensitivity_control: {pass_sens}",
        f"- orthogonality_proxy_corr: {pass_corr}",
        f"- spread_control: {pass_spread}",
        f"- nonboundary_omega: {pass_nonboundary}",
        f"- fit_quality: {pass_fit}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1756] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1756] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
