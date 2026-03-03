#!/usr/bin/env python3
"""
QW-1802: Decoherence form comparison (multiplicative vs additive).

Follow-up to QW-1801 (best damping kind: quadratic).
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
from scipy.special import logsumexp
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1802_decoherence_form_comparison.json"
OUT_MD = ROOT / "RAPORT_QW1802_DECOHERENCE_FORM_COMPARISON.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def damping_quadratic(theta_deg: np.ndarray, lam: float) -> np.ndarray:
    x = (1.0 - np.cos(np.radians(theta_deg))) / 2.0
    return np.exp(-lam * (x * x))


def additive_feature(theta_deg: np.ndarray, lam: float) -> np.ndarray:
    d = damping_quadratic(theta_deg, lam)
    g = 1.0 - d
    return g - float(np.mean(g))


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


def evidence_m6_mul(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    lam = rng.uniform(0.0, 6.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    cache: Dict[float, np.ndarray] = {}
    ll_rows = []
    for a, c, qq, ll in zip(A, C, q, lam):
        key = float(np.round(ll, 3))
        if key not in cache:
            cache[key] = damping_quadratic(theta, key)
        d = cache[key]
        ll_rows.append(loglike(H, sigma, a * (hd0 ** qq) * d + c))
    ll = np.array(ll_rows, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_m6_add(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.2, 1.2, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    lam = rng.uniform(0.0, 6.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    cache: Dict[float, np.ndarray] = {}
    ll_rows = []
    for a, b, c, qq, ll in zip(A, B, C, q, lam):
        key = float(np.round(ll, 3))
        if key not in cache:
            cache[key] = additive_feature(theta, key)
        g = cache[key]
        ll_rows.append(loglike(H, sigma, a * (hd0 ** qq) + b * g + c))
    ll = np.array(ll_rows, dtype=float)
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


def fit_best_m6_mul(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    lam = rng.uniform(0.0, 6.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    cache: Dict[float, np.ndarray] = {}
    best_ll = -np.inf
    best = None
    for a, c, qq, ll in zip(A, C, q, lam):
        key = float(np.round(ll, 3))
        if key not in cache:
            cache[key] = damping_quadratic(theta, key)
        d = cache[key]
        m = a * (hd0 ** qq) * d + c
        lk = loglike(H, sigma, m)
        if lk > best_ll:
            best_ll = lk
            best = {"A": float(a), "C": float(c), "q": float(qq), "lambda": float(ll), "sigma": sigma}
    return best


def fit_best_m6_add(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.2, 1.2, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    lam = rng.uniform(0.0, 6.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    cache: Dict[float, np.ndarray] = {}
    best_ll = -np.inf
    best = None
    for a, b, c, qq, ll in zip(A, B, C, q, lam):
        key = float(np.round(ll, 3))
        if key not in cache:
            cache[key] = additive_feature(theta, key)
        g = cache[key]
        m = a * (hd0 ** qq) + b * g + c
        lk = loglike(H, sigma, m)
        if lk > best_ll:
            best_ll = lk
            best = {"A": float(a), "B": float(b), "C": float(c), "q": float(qq), "lambda": float(ll), "sigma": sigma}
    return best


def residual_stats(resid: np.ndarray, theta: np.ndarray) -> Dict[str, float]:
    rho, _ = spearmanr(theta, resid)
    rho = float(rho)
    bins = np.linspace(0.0, 180.0, 7)
    means = []
    for i in range(len(bins) - 1):
        m = (theta >= bins[i]) & (theta < bins[i + 1] if i < len(bins) - 2 else theta <= bins[i + 1])
        if np.sum(m) == 0:
            continue
        means.append(float(np.mean(resid[m])))
    max_bin = float(max(abs(v) for v in means)) if means else 0.0
    return {"rho_theta_resid": rho, "max_abs_bin_mean": max_bin}


def main() -> None:
    helper = load_module(ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py", "qw1787_helper")
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
        stab = helper.split_half_stability(x, y)
        if len(x) >= n_match_min and stab <= stability_max:
            rows.append({"theta_deg": float(sep), "hxy": float(hxy)})

    if len(rows) < 85:
        raise RuntimeError("Operational cohort too small.")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    hd0 = np.clip(helper.hellings_downs(theta_all), 1e-9, None)

    # Full evidence
    z0 = evidence_flat(H_all, n_mc=8500, seed=18021)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=18022)
    z6m = evidence_m6_mul(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=17000, seed=18023)
    z6a = evidence_m6_add(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=19000, seed=18024)

    full_m2 = float(z2 - z0)
    full_m6m = float(z6m - z0)
    full_m6a = float(z6a - z0)
    full_delta_m6m = float(full_m6m - full_m2)
    full_delta_m6a = float(full_m6a - full_m2)
    full_delta_add_minus_mul = float(full_m6a - full_m6m)

    # Residual diagnostics
    fit2 = fit_best_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=60000, seed=18025)
    fit6m = fit_best_m6_mul(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=85000, seed=18026)
    fit6a = fit_best_m6_add(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=90000, seed=18027)

    pred2 = fit2["A"] * (hd0 ** fit2["q"]) + fit2["C"]
    d6 = damping_quadratic(theta_all, fit6m["lambda"])
    pred6m = fit6m["A"] * (hd0 ** fit6m["q"]) * d6 + fit6m["C"]
    g6 = additive_feature(theta_all, fit6a["lambda"])
    pred6a = fit6a["A"] * (hd0 ** fit6a["q"]) + fit6a["B"] * g6 + fit6a["C"]

    rs2 = residual_stats(H_all - pred2, theta_all)
    rsm = residual_stats(H_all - pred6m, theta_all)
    rsa = residual_stats(H_all - pred6a, theta_all)

    # Replications
    rng = np.random.default_rng(18028)
    reps = []
    for i in range(12):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 80:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        b0 = evidence_flat(hh, n_mc=4000, seed=18030 + 20 * i + 1)
        b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=18030 + 20 * i + 2)
        b6m = evidence_m6_mul(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=7800, seed=18030 + 20 * i + 3)
        b6a = evidence_m6_add(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=8600, seed=18030 + 20 * i + 4)
        l2 = float(b2 - b0)
        l6m = float(b6m - b0)
        l6a = float(b6a - b0)
        reps.append(
            {
                "rep": i,
                "logB_m2_vs_flat": l2,
                "logB_m6_mul_vs_flat": l6m,
                "logB_m6_add_vs_flat": l6a,
                "delta_mul_vs_m2": float(l6m - l2),
                "delta_add_vs_m2": float(l6a - l2),
                "delta_add_vs_mul": float(l6a - l6m),
            }
        )

    arr_add = np.array([r["logB_m6_add_vs_flat"] for r in reps], dtype=float)
    arr_mul = np.array([r["logB_m6_mul_vs_flat"] for r in reps], dtype=float)
    arr_dm = np.array([r["delta_mul_vs_m2"] for r in reps], dtype=float)
    arr_da = np.array([r["delta_add_vs_m2"] for r in reps], dtype=float)
    arr_am = np.array([r["delta_add_vs_mul"] for r in reps], dtype=float)

    summary = {
        "n_pairs": len(rows),
        "fraction": frac,
        "q_width": q_width,
        "full_logB_m2_vs_flat": full_m2,
        "full_logB_m6_mul_vs_flat": full_m6m,
        "full_logB_m6_add_vs_flat": full_m6a,
        "full_delta_mul_vs_m2": full_delta_m6m,
        "full_delta_add_vs_m2": full_delta_m6a,
        "full_delta_add_vs_mul": full_delta_add_minus_mul,
        "replications": len(reps),
        "prob_add_gt_mul": float(np.mean(arr_am > 0.0)),
        "prob_add_gt_m2": float(np.mean(arr_da > 0.0)),
        "prob_add_gt_flat": float(np.mean(arr_add > 0.0)),
        "std_delta_add_vs_m2": float(np.std(arr_da)),
        "prob_mul_gt_m2": float(np.mean(arr_dm > 0.0)),
        "prob_mul_gt_flat": float(np.mean(arr_mul > 0.0)),
        "std_delta_mul_vs_m2": float(np.std(arr_dm)),
        "resid_m2_rho_theta": rs2["rho_theta_resid"],
        "resid_mul_rho_theta": rsm["rho_theta_resid"],
        "resid_add_rho_theta": rsa["rho_theta_resid"],
        "resid_m2_max_abs_bin_mean": rs2["max_abs_bin_mean"],
        "resid_mul_max_abs_bin_mean": rsm["max_abs_bin_mean"],
        "resid_add_max_abs_bin_mean": rsa["max_abs_bin_mean"],
        "resid_add_rho_improvement_abs": float(abs(rs2["rho_theta_resid"]) - abs(rsa["rho_theta_resid"])),
        "resid_add_bin_improvement": float(rs2["max_abs_bin_mean"] - rsa["max_abs_bin_mean"]),
    }

    pass_add_gain = summary["full_delta_add_vs_m2"] > 0.0 and summary["prob_add_gt_m2"] >= 0.75 and summary["prob_add_gt_flat"] >= 0.95
    pass_add_vs_mul = summary["full_delta_add_vs_mul"] > 0.0 and summary["prob_add_gt_mul"] >= 0.65
    pass_add_disp = summary["std_delta_add_vs_m2"] <= 0.25
    pass_add_resid = summary["resid_add_rho_improvement_abs"] > 0.03 and summary["resid_add_bin_improvement"] > 0.01

    if pass_add_gain and pass_add_vs_mul and pass_add_disp and pass_add_resid:
        verdict = "DECOHERENCE_ADDITIVE_FORM_SUPPORTED"
    elif pass_add_gain and pass_add_disp:
        verdict = "DECOHERENCE_ADDITIVE_FORM_PARTIAL"
    else:
        verdict = "DECOHERENCE_ADDITIVE_FORM_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "fit_params": {"m2": fit2, "m6_mul": fit6m, "m6_add": fit6a},
        "pass_flags": {
            "additive_gain_over_m2": bool(pass_add_gain),
            "additive_better_than_multiplicative": bool(pass_add_vs_mul),
            "additive_dispersion_control": bool(pass_add_disp),
            "additive_residual_flattening": bool(pass_add_resid),
        },
        "verdict": verdict,
        "replications": reps,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1802: DECOHERENCE FORM COMPARISON",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Full logB M2/mul/add: {full_m2:.4f} / {full_m6m:.4f} / {full_m6a:.4f}",
        f"- Full delta add-vs-M2: {full_delta_m6a:.4f}",
        f"- Full delta add-vs-mul: {full_delta_add_minus_mul:.4f}",
        f"- P(add>mul): {summary['prob_add_gt_mul']:.3f}",
        f"- P(add>M2): {summary['prob_add_gt_m2']:.3f}",
        f"- P(add>flat): {summary['prob_add_gt_flat']:.3f}",
        f"- Std delta add-vs-M2: {summary['std_delta_add_vs_m2']:.3f}",
        f"- Residual add improvements: d|rho|={summary['resid_add_rho_improvement_abs']:.4f}, d(bin)={summary['resid_add_bin_improvement']:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- additive_gain_over_m2: {pass_add_gain}",
        f"- additive_better_than_multiplicative: {pass_add_vs_mul}",
        f"- additive_dispersion_control: {pass_add_disp}",
        f"- additive_residual_flattening: {pass_add_resid}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1802] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1802] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
