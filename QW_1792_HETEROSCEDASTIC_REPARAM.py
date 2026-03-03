#!/usr/bin/env python3
"""
QW-1792: Heteroscedastic-noise reparam test.

Rationale:
- keep the same signal family (flat / legacy / reparam),
- upgrade only the measurement layer to heteroscedastic Gaussian noise
  using pair-quality proxies (n_match, split-half stability).
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


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1792_heteroscedastic_reparam.json"
OUT_MD = ROOT / "RAPORT_QW1792_HETEROSCEDASTIC_REPARAM.md"


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
    # Higher proxy -> lower quality -> larger effective uncertainty.
    inv_sqrt_n = 1.0 / np.sqrt(np.maximum(n_match, 1.0))
    a = minmax(inv_sqrt_n)
    b = minmax(stability)
    return 0.5 * (a + b)


def loglike_full(H: np.ndarray, sigma_i: np.ndarray, model: np.ndarray) -> float:
    var = sigma_i * sigma_i
    resid2 = (H - model) ** 2
    return float(-0.5 * np.sum(resid2 / var + np.log(2.0 * np.pi * var)))


def sigma_from_eta(base_sigma: float, proxy: np.ndarray, eta: float) -> np.ndarray:
    # Keep average scale close to base_sigma while allowing quality-weighted spread.
    centered = proxy - float(np.mean(proxy))
    return base_sigma * np.exp(eta * centered)


def evidence_homo_flat(H: np.ndarray, base_sigma: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 2.0, n_mc)
    sig = np.full_like(H, base_sigma, dtype=float)
    ll = np.array([loglike_full(H, sig, np.full_like(H, c)) for c in C], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_homo_legacy(helper, theta: np.ndarray, H: np.ndarray, base_sigma: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    hd = helper.hellings_downs(theta)
    sig = np.full_like(H, base_sigma, dtype=float)
    ll = np.array([loglike_full(H, sig, a * hd + c) for a, c in zip(A, C)], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_homo_reparam(helper, theta: np.ndarray, H: np.ndarray, base_sigma: float, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    sig = np.full_like(H, base_sigma, dtype=float)
    ll = np.array([loglike_full(H, sig, a * (hd0 ** qq) + c) for a, c, qq in zip(A, C, q)], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_hetero_flat(H: np.ndarray, base_sigma: float, proxy: np.ndarray, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 2.0, n_mc)
    eta = rng.uniform(0.0, 1.0, n_mc)
    ll = np.array(
        [loglike_full(H, sigma_from_eta(base_sigma, proxy, e), np.full_like(H, c)) for c, e in zip(C, eta)],
        dtype=float,
    )
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_hetero_legacy(helper, theta: np.ndarray, H: np.ndarray, base_sigma: float, proxy: np.ndarray, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    eta = rng.uniform(0.0, 1.0, n_mc)
    hd = helper.hellings_downs(theta)
    ll = np.array(
        [loglike_full(H, sigma_from_eta(base_sigma, proxy, e), a * hd + c) for a, c, e in zip(A, C, eta)],
        dtype=float,
    )
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_hetero_reparam(helper, theta: np.ndarray, H: np.ndarray, base_sigma: float, proxy: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    eta = rng.uniform(0.0, 1.0, n_mc)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array(
        [
            loglike_full(H, sigma_from_eta(base_sigma, proxy, e), a * (hd0 ** qq) + c)
            for a, c, qq, e in zip(A, C, q, eta)
        ],
        dtype=float,
    )
    return float(logsumexp(ll) - np.log(len(ll)))


def run_family(helper, theta: np.ndarray, H: np.ndarray, proxy: np.ndarray, q_center: float, q_width: float, mode: str, seed_base: int) -> Dict[str, float]:
    base_sigma = max(float(np.std(H)), 1e-6)
    if mode == "homo":
        z0 = evidence_homo_flat(H, base_sigma, n_mc=9000, seed=seed_base + 1)
        z1 = evidence_homo_legacy(helper, theta, H, base_sigma, n_mc=12000, seed=seed_base + 2)
        z2 = evidence_homo_reparam(helper, theta, H, base_sigma, q_center=q_center, q_width=q_width, n_mc=16000, seed=seed_base + 3)
    else:
        z0 = evidence_hetero_flat(H, base_sigma, proxy, n_mc=10000, seed=seed_base + 1)
        z1 = evidence_hetero_legacy(helper, theta, H, base_sigma, proxy, n_mc=13000, seed=seed_base + 2)
        z2 = evidence_hetero_reparam(helper, theta, H, base_sigma, proxy, q_center=q_center, q_width=q_width, n_mc=17000, seed=seed_base + 3)
    return {
        "logB_legacy_vs_flat": float(z1 - z0),
        "logB_reparam_vs_flat": float(z2 - z0),
        "delta_logB_reparam_minus_legacy": float((z2 - z0) - (z1 - z0)),
    }


def main() -> None:
    helper = load_helper()

    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    q_center = float(d1773["projection"]["p"])
    q_width = 0.20
    frac = 0.95
    n_rep = 14

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
        if len(x) >= 120 and stab <= 0.65:
            rows.append(
                {
                    "theta_deg": float(sep),
                    "hxy": float(hxy),
                    "n_match": float(len(x)),
                    "stability": float(stab),
                }
            )

    if len(rows) < 85:
        raise RuntimeError("Base operational cohort too small.")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    n_match_all = np.array([r["n_match"] for r in rows], dtype=float)
    stab_all = np.array([r["stability"] for r in rows], dtype=float)
    proxy_all = quality_proxy(n_match_all, stab_all)

    full_homo = run_family(helper, theta_all, H_all, proxy_all, q_center=q_center, q_width=q_width, mode="homo", seed_base=14000)
    full_hetero = run_family(helper, theta_all, H_all, proxy_all, q_center=q_center, q_width=q_width, mode="hetero", seed_base=15000)

    rng = np.random.default_rng(1792)
    rep_rows = []
    for i in range(n_rep):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 80:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        proxy = proxy_all[idx]
        r_homo = run_family(helper, th, hh, proxy, q_center=q_center, q_width=q_width, mode="homo", seed_base=16000 + 30 * i)
        r_hetero = run_family(helper, th, hh, proxy, q_center=q_center, q_width=q_width, mode="hetero", seed_base=17000 + 30 * i)
        rep_rows.append(
            {
                "rep": i,
                "n_pairs": int(len(idx)),
                "homo_logB_reparam_vs_flat": r_homo["logB_reparam_vs_flat"],
                "hetero_logB_reparam_vs_flat": r_hetero["logB_reparam_vs_flat"],
                "delta_hetero_minus_homo_reparam_vs_flat": float(r_hetero["logB_reparam_vs_flat"] - r_homo["logB_reparam_vs_flat"]),
                "homo_delta_reparam_minus_legacy": r_homo["delta_logB_reparam_minus_legacy"],
                "hetero_delta_reparam_minus_legacy": r_hetero["delta_logB_reparam_minus_legacy"],
            }
        )

    arr_hetero = np.array([r["hetero_logB_reparam_vs_flat"] for r in rep_rows], dtype=float)
    arr_gain = np.array([r["delta_hetero_minus_homo_reparam_vs_flat"] for r in rep_rows], dtype=float)
    arr_delta_rl = np.array([r["hetero_delta_reparam_minus_legacy"] for r in rep_rows], dtype=float)

    summary = {
        "n_pairs_operational": len(rows),
        "q_width": q_width,
        "fraction": frac,
        "replications": len(rep_rows),
        "full_homo_logB_reparam_vs_flat": full_homo["logB_reparam_vs_flat"],
        "full_hetero_logB_reparam_vs_flat": full_hetero["logB_reparam_vs_flat"],
        "full_gain_hetero_minus_homo_reparam_vs_flat": float(full_hetero["logB_reparam_vs_flat"] - full_homo["logB_reparam_vs_flat"]),
        "full_homo_delta_reparam_minus_legacy": full_homo["delta_logB_reparam_minus_legacy"],
        "full_hetero_delta_reparam_minus_legacy": full_hetero["delta_logB_reparam_minus_legacy"],
        "prob_hetero_positive": float(np.mean(arr_hetero > 0.0)),
        "prob_hetero_better_than_homo": float(np.mean(arr_gain > 0.0)),
        "prob_hetero_delta_vs_legacy_positive": float(np.mean(arr_delta_rl > 0.0)),
        "median_gain_hetero_minus_homo": float(np.median(arr_gain)),
        "std_gain_hetero_minus_homo": float(np.std(arr_gain)),
    }

    pass_full_positive = summary["full_hetero_logB_reparam_vs_flat"] > 0.0
    pass_full_gain = summary["full_gain_hetero_minus_homo_reparam_vs_flat"] > 0.0
    pass_rep_positive = summary["prob_hetero_positive"] >= 0.95
    pass_rep_gain = summary["prob_hetero_better_than_homo"] >= 0.70
    pass_rep_delta = summary["prob_hetero_delta_vs_legacy_positive"] >= 0.95
    pass_disp = summary["std_gain_hetero_minus_homo"] <= 0.16

    if pass_full_positive and pass_full_gain and pass_rep_positive and pass_rep_gain and pass_rep_delta and pass_disp:
        verdict = "HETEROSCEDASTIC_REPARAM_SUPPORTED"
    elif pass_full_positive and pass_rep_positive and pass_rep_delta:
        verdict = "HETEROSCEDASTIC_REPARAM_PARTIAL"
    else:
        verdict = "HETEROSCEDASTIC_REPARAM_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "pass_flags": {
            "full_positive": bool(pass_full_positive),
            "full_gain_over_homo": bool(pass_full_gain),
            "rep_positive": bool(pass_rep_positive),
            "rep_gain_over_homo": bool(pass_rep_gain),
            "rep_delta_vs_legacy_positive": bool(pass_rep_delta),
            "dispersion_control": bool(pass_disp),
        },
        "verdict": verdict,
        "full_homoscedastic": full_homo,
        "full_heteroscedastic": full_hetero,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1792: HETEROSCEDASTIC REPARAM",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Operational pairs: {len(rows)}",
        f"- Full homo logB(reparam vs flat): {summary['full_homo_logB_reparam_vs_flat']:.4f}",
        f"- Full hetero logB(reparam vs flat): {summary['full_hetero_logB_reparam_vs_flat']:.4f}",
        f"- Full gain hetero-homo: {summary['full_gain_hetero_minus_homo_reparam_vs_flat']:.4f}",
        f"- P(hetero>0): {summary['prob_hetero_positive']:.3f}",
        f"- P(hetero>homo): {summary['prob_hetero_better_than_homo']:.3f}",
        f"- P(hetero delta vs legacy > 0): {summary['prob_hetero_delta_vs_legacy_positive']:.3f}",
        f"- Std gain hetero-homo: {summary['std_gain_hetero_minus_homo']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_positive: {pass_full_positive}",
        f"- full_gain_over_homo: {pass_full_gain}",
        f"- rep_positive: {pass_rep_positive}",
        f"- rep_gain_over_homo: {pass_rep_gain}",
        f"- rep_delta_vs_legacy_positive: {pass_rep_delta}",
        f"- dispersion_control: {pass_disp}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1792] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1792] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
