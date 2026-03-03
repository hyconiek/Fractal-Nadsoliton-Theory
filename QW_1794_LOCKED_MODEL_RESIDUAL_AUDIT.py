#!/usr/bin/env python3
"""
QW-1794: Residual-structure audit for locked operational model.

Locked model from QW-1793:
    H(theta) = A * HD(theta)^q + C
with homoscedastic Gaussian likelihood.
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1794_locked_model_residual_audit.json"
OUT_MD = ROOT / "RAPORT_QW1794_LOCKED_MODEL_RESIDUAL_AUDIT.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def minmax(v: np.ndarray) -> np.ndarray:
    lo = float(np.min(v))
    hi = float(np.max(v))
    if hi <= lo:
        return np.zeros_like(v)
    return (v - lo) / (hi - lo)


def quality_proxy(n_match: np.ndarray, stability: np.ndarray) -> np.ndarray:
    inv_sqrt_n = 1.0 / np.sqrt(np.maximum(n_match, 1.0))
    return 0.5 * (minmax(inv_sqrt_n) + minmax(stability))


def best_fit_m2(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_samples: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_samples)
    C = rng.uniform(-1.0, 2.0, n_samples)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_samples), 0.8, 2.4)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    sigma = max(float(np.std(H)), 1e-6)

    best = None
    best_ll = -np.inf
    for a, c, qq in zip(A, C, q):
        m = a * (hd0 ** qq) + c
        z = (H - m) / sigma
        ll = float(-0.5 * np.sum(z * z))
        if ll > best_ll:
            best_ll = ll
            best = {"A": float(a), "C": float(c), "q": float(qq), "sigma": sigma, "loglike": best_ll}
    return best


def permutation_pvalue(x: np.ndarray, y: np.ndarray, obs: float, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        yp = rng.permutation(y)
        r, _ = spearmanr(x, yp)
        if abs(float(r)) >= abs(obs):
            count += 1
    return float((count + 1) / (n_perm + 1))


def main() -> None:
    helper = load_helper()
    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))

    q_center = float(d1773["projection"]["p"])
    q_width = float(d1793["operational_protocol"]["q_width"])
    cohort = d1793["operational_protocol"]["cohort"]
    n_match_min = int(cohort["n_match_min"])
    stability_max = float(cohort["stability_max"])

    residuals = helper.load_residuals(ROOT / "nano15/residuals/NANOGrav15yr_PulsarTiming_v2.1.0/residuals", max_psr=34)
    positions = helper.load_positions(ROOT / "nano15/parfiles")

    rows: List[Dict[str, float]] = []
    psr_list = list(residuals.keys())
    for p1, p2 in combinations(psr_list, 2):
        sep = helper.angular_sep(p1, p2, positions)
        if sep is None:
            continue
        x, y = helper.match_epochs(residuals[p1], residuals[p2], tol_days=30.0)
        if x is None:
            continue
        hxy = helper.cross_dfa(x, y, min_scale=15)
        if not np.isfinite(hxy):
            continue
        stab = helper.split_half_stability(x, y)
        if len(x) >= n_match_min and stab <= stability_max:
            rows.append(
                {
                    "theta_deg": float(sep),
                    "hxy": float(hxy),
                    "n_match": float(len(x)),
                    "stability": float(stab),
                }
            )

    if len(rows) < 85:
        raise RuntimeError("Operational cohort too small.")

    theta = np.array([r["theta_deg"] for r in rows], dtype=float)
    H = np.array([r["hxy"] for r in rows], dtype=float)
    n_match = np.array([r["n_match"] for r in rows], dtype=float)
    stability = np.array([r["stability"] for r in rows], dtype=float)
    qproxy = quality_proxy(n_match, stability)

    fit = best_fit_m2(helper, theta, H, q_center=q_center, q_width=q_width, n_samples=120000, seed=17940)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    pred = fit["A"] * (hd0 ** fit["q"]) + fit["C"]
    resid = H - pred

    rho_theta, _ = spearmanr(theta, resid)
    rho_quality, _ = spearmanr(qproxy, resid)
    rho_theta = float(rho_theta)
    rho_quality = float(rho_quality)
    p_theta = permutation_pvalue(theta, resid, rho_theta, n_perm=1500, seed=17941)
    p_quality = permutation_pvalue(qproxy, resid, rho_quality, n_perm=1500, seed=17942)

    # Angle-binned residual means.
    bins = np.linspace(0.0, 180.0, 7)
    bin_stats = []
    for i in range(len(bins) - 1):
        m = (theta >= bins[i]) & (theta < bins[i + 1] if i < len(bins) - 2 else theta <= bins[i + 1])
        if np.sum(m) == 0:
            continue
        rr = resid[m]
        bin_stats.append(
            {
                "bin_deg": [float(bins[i]), float(bins[i + 1])],
                "n": int(np.sum(m)),
                "mean_residual": float(np.mean(rr)),
                "std_residual": float(np.std(rr)),
            }
        )
    max_abs_bin_mean = float(max(abs(b["mean_residual"]) for b in bin_stats))

    # Standardized residual diagnostics.
    sigma = float(fit["sigma"])
    z = resid / max(sigma, 1e-9)
    z_mean = float(np.mean(z))
    z_std = float(np.std(z))
    z_q95 = float(np.quantile(np.abs(z), 0.95))

    pass_theta = abs(rho_theta) <= 0.15 and p_theta >= 0.05
    pass_quality = abs(rho_quality) <= 0.15 and p_quality >= 0.05
    pass_bins = max_abs_bin_mean <= 0.08
    pass_z = abs(z_mean) <= 0.05 and 0.85 <= z_std <= 1.20 and z_q95 <= 2.8

    if pass_theta and pass_quality and pass_bins and pass_z:
        verdict = "LOCKED_MODEL_RESIDUALS_NO_STRONG_STRUCTURE"
    elif pass_theta and pass_quality:
        verdict = "LOCKED_MODEL_RESIDUALS_PARTIAL_STRUCTURE"
    else:
        verdict = "LOCKED_MODEL_RESIDUALS_STRUCTURED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cohort_size": len(rows),
        "fit": fit,
        "diagnostics": {
            "rho_theta_residual": rho_theta,
            "rho_quality_residual": rho_quality,
            "perm_p_theta": p_theta,
            "perm_p_quality": p_quality,
            "max_abs_bin_mean_residual": max_abs_bin_mean,
            "z_mean": z_mean,
            "z_std": z_std,
            "z_abs_q95": z_q95,
        },
        "pass_flags": {
            "theta_independence": bool(pass_theta),
            "quality_independence": bool(pass_quality),
            "angle_bin_flatness": bool(pass_bins),
            "z_residual_diagnostics": bool(pass_z),
        },
        "verdict": verdict,
        "bin_stats": bin_stats,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1794: LOCKED MODEL RESIDUAL AUDIT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Cohort size: {len(rows)}",
        f"- Best-fit params: A={fit['A']:.4f}, q={fit['q']:.4f}, C={fit['C']:.4f}, sigma={fit['sigma']:.4f}",
        f"- rho(theta,resid)={rho_theta:.4f}, p_perm={p_theta:.4f}",
        f"- rho(quality,resid)={rho_quality:.4f}, p_perm={p_quality:.4f}",
        f"- max |bin mean residual|={max_abs_bin_mean:.4f}",
        f"- z diagnostics: mean={z_mean:.4f}, std={z_std:.4f}, q95|z|={z_q95:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- theta_independence: {pass_theta}",
        f"- quality_independence: {pass_quality}",
        f"- angle_bin_flatness: {pass_bins}",
        f"- z_residual_diagnostics: {pass_z}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1794] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1794] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
