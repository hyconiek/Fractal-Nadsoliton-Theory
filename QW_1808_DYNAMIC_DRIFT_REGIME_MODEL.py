#!/usr/bin/env python3
"""
QW-1808: Dynamic drift-regime model (phase-1 dynamic program).

Dynamic latent proxy:
    drift = hxy(second_half) - hxy(first_half)

Model:
    M2D(theta, drift) = A * HD(theta)^q + B * drift + C

This is the first dynamic extension after QW-1807 transition gate.
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.special import logsumexp
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1808_dynamic_drift_regime_model.json"
OUT_MD = ROOT / "RAPORT_QW1808_DYNAMIC_DRIFT_REGIME_MODEL.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def split_half_hxy(helper, x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    n = min(len(x), len(y))
    if n < 140:
        return float("nan"), float("nan")
    h = n // 2
    h1 = helper.cross_dfa(x[:h], y[:h], min_scale=12)
    h2 = helper.cross_dfa(x[h:], y[h:], min_scale=12)
    return float(h1), float(h2)


def loglike(H: np.ndarray, sigma: float, model: np.ndarray) -> float:
    z = (H - model) / sigma
    return float(-0.5 * np.sum(z * z))


def evidence_flat(H: np.ndarray, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 2.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    ll = np.array([loglike(H, sigma, np.full_like(H, c)) for c in C], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_m2(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array([loglike(H, sigma, a * (hd0 ** qq) + c) for a, c, qq in zip(A, C, q)], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_m2d(helper, theta: np.ndarray, H: np.ndarray, D: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.2, 1.2, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array([loglike(H, sigma, a * (hd0 ** qq) + b * D + c) for a, b, c, qq in zip(A, B, C, q)], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def fit_best_m2(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    best_ll = -np.inf
    best = None
    for a, c, qq in zip(A, C, q):
        m = a * (hd0 ** qq) + c
        lk = loglike(H, sigma, m)
        if lk > best_ll:
            best_ll = lk
            best = {"A": float(a), "C": float(c), "q": float(qq), "sigma": sigma}
    return best


def fit_best_m2d(helper, theta: np.ndarray, H: np.ndarray, D: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.2, 1.2, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    best_ll = -np.inf
    best = None
    for a, b, c, qq in zip(A, B, C, q):
        m = a * (hd0 ** qq) + b * D + c
        lk = loglike(H, sigma, m)
        if lk > best_ll:
            best_ll = lk
            best = {"A": float(a), "B": float(b), "C": float(c), "q": float(qq), "sigma": sigma}
    return best


def main() -> None:
    helper = load_helper()

    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))

    q_center = float(d1773["projection"]["p"])
    q_width = float(d1793["operational_protocol"]["q_width"])
    frac = float(d1793["operational_protocol"]["fraction"])
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
        h1, h2 = split_half_hxy(helper, x, y)
        if not np.isfinite(h1) or not np.isfinite(h2):
            continue
        drift = float(h2 - h1)
        stability = float(abs(drift))
        if len(x) >= n_match_min and stability <= stability_max:
            rows.append(
                {
                    "theta_deg": float(sep),
                    "hxy": float(hxy),
                    "drift": drift,
                }
            )

    if len(rows) < 85:
        raise RuntimeError("Operational cohort too small.")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    D_all = np.array([r["drift"] for r in rows], dtype=float)
    D_all = D_all - float(np.mean(D_all))
    hd0 = np.clip(helper.hellings_downs(theta_all), 1e-9, None)

    z0 = evidence_flat(H_all, n_mc=8500, seed=18801)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=18802)
    z2d = evidence_m2d(helper, theta_all, H_all, D_all, q_center=q_center, q_width=q_width, n_mc=17000, seed=18803)
    full_m2 = float(z2 - z0)
    full_m2d = float(z2d - z0)
    full_delta = float(full_m2d - full_m2)

    fit2 = fit_best_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=65000, seed=18804)
    fit2d = fit_best_m2d(helper, theta_all, H_all, D_all, q_center=q_center, q_width=q_width, n_mc=85000, seed=18805)
    pred2 = fit2["A"] * (hd0 ** fit2["q"]) + fit2["C"]
    pred2d = fit2d["A"] * (hd0 ** fit2d["q"]) + fit2d["B"] * D_all + fit2d["C"]
    resid2 = H_all - pred2
    resid2d = H_all - pred2d
    rho_d_resid2 = float(spearmanr(D_all, resid2)[0])
    rho_d_resid2d = float(spearmanr(D_all, resid2d)[0])
    drift_improvement = float(abs(rho_d_resid2) - abs(rho_d_resid2d))

    rng = np.random.default_rng(18806)
    rep_rows = []
    for i in range(12):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 80:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        dd = D_all[idx]
        b0 = evidence_flat(hh, n_mc=3900, seed=18810 + 20 * i + 1)
        b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=18810 + 20 * i + 2)
        b2d = evidence_m2d(helper, th, hh, dd, q_center=q_center, q_width=q_width, n_mc=8200, seed=18810 + 20 * i + 3)
        l2 = float(b2 - b0)
        l2d = float(b2d - b0)
        rep_rows.append({"rep": i, "logB_m2_vs_flat": l2, "logB_m2d_vs_flat": l2d, "delta_m2d_vs_m2": float(l2d - l2)})

    arr_d = np.array([r["delta_m2d_vs_m2"] for r in rep_rows], dtype=float)
    arr_m = np.array([r["logB_m2d_vs_flat"] for r in rep_rows], dtype=float)

    summary = {
        "n_pairs": len(rows),
        "fraction": frac,
        "q_width": q_width,
        "full_logB_m2_vs_flat": full_m2,
        "full_logB_m2d_vs_flat": full_m2d,
        "full_delta_m2d_vs_m2": full_delta,
        "replications": len(rep_rows),
        "prob_m2d_gt_m2": float(np.mean(arr_d > 0.0)),
        "prob_m2d_gt_flat": float(np.mean(arr_m > 0.0)),
        "median_delta_m2d_vs_m2": float(np.median(arr_d)),
        "std_delta_m2d_vs_m2": float(np.std(arr_d)),
        "rho_drift_resid_m2": rho_d_resid2,
        "rho_drift_resid_m2d": rho_d_resid2d,
        "drift_residual_improvement_abs": drift_improvement,
    }

    pass_full = summary["full_delta_m2d_vs_m2"] > 0.0
    pass_rep = summary["prob_m2d_gt_m2"] >= 0.80 and summary["prob_m2d_gt_flat"] >= 0.95
    pass_disp = summary["std_delta_m2d_vs_m2"] <= 0.30
    pass_dyn = summary["drift_residual_improvement_abs"] > 0.05

    if pass_full and pass_rep and pass_disp and pass_dyn:
        verdict = "DYNAMIC_DRIFT_REGIME_SUPPORTED"
    elif pass_full and pass_rep and pass_dyn:
        verdict = "DYNAMIC_DRIFT_REGIME_PARTIAL"
    else:
        verdict = "DYNAMIC_DRIFT_REGIME_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "fit_params": {"m2": fit2, "m2d": fit2d},
        "pass_flags": {
            "full_gain": bool(pass_full),
            "replication_gain": bool(pass_rep),
            "dispersion_control": bool(pass_disp),
            "dynamic_residual_reduction": bool(pass_dyn),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1808: DYNAMIC DRIFT REGIME MODEL",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Full logB M2/M2D: {full_m2:.4f} / {full_m2d:.4f}",
        f"- Full delta M2D-M2: {full_delta:.4f}",
        f"- P(M2D>M2): {summary['prob_m2d_gt_m2']:.3f}",
        f"- P(M2D>flat): {summary['prob_m2d_gt_flat']:.3f}",
        f"- Std delta M2D-M2: {summary['std_delta_m2d_vs_m2']:.3f}",
        f"- Dynamic residual improvement |rho|: {drift_improvement:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_gain: {pass_full}",
        f"- replication_gain: {pass_rep}",
        f"- dispersion_control: {pass_disp}",
        f"- dynamic_residual_reduction: {pass_dyn}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1808] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1808] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
