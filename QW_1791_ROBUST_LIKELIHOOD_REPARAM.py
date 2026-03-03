#!/usr/bin/env python3
"""
QW-1791: Robust-likelihood test for reparametrized PTA model.

Compares:
- Gaussian likelihood (baseline)
- Student-t likelihood (robust, nu fixed)
for flat / legacy / reparam models on the same cohort and protocol.
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
OUT_JSON = ROOT / "report_qw1791_robust_likelihood_reparam.json"
OUT_MD = ROOT / "RAPORT_QW1791_ROBUST_LIKELIHOOD_REPARAM.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def loglike_gaussian(H: np.ndarray, sigma: float, model: np.ndarray) -> float:
    z = (H - model) / sigma
    return float(-0.5 * np.sum(z * z))


def loglike_student_t(H: np.ndarray, sigma: float, model: np.ndarray, nu: float) -> float:
    z = (H - model) / sigma
    return float(-0.5 * (nu + 1.0) * np.sum(np.log1p((z * z) / nu)))


def evidence_flat(H: np.ndarray, sigma: float, n_mc: int, seed: int, like_mode: str, nu: float) -> float:
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 2.0, n_mc)
    ll = []
    for c in C:
        m = np.full_like(H, c)
        if like_mode == "gauss":
            ll.append(loglike_gaussian(H, sigma, m))
        else:
            ll.append(loglike_student_t(H, sigma, m, nu))
    ll = np.array(ll, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_legacy(helper, theta: np.ndarray, H: np.ndarray, sigma: float, n_mc: int, seed: int, like_mode: str, nu: float) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    hd = helper.hellings_downs(theta)
    ll = []
    for a, c in zip(A, C):
        m = a * hd + c
        if like_mode == "gauss":
            ll.append(loglike_gaussian(H, sigma, m))
        else:
            ll.append(loglike_student_t(H, sigma, m, nu))
    ll = np.array(ll, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_reparam(helper, theta: np.ndarray, H: np.ndarray, sigma: float, q_center: float, q_width: float, n_mc: int, seed: int, like_mode: str, nu: float) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = []
    for a, c, qq in zip(A, C, q):
        m = a * (hd0 ** qq) + c
        if like_mode == "gauss":
            ll.append(loglike_gaussian(H, sigma, m))
        else:
            ll.append(loglike_student_t(H, sigma, m, nu))
    ll = np.array(ll, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def run_family(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, like_mode: str, nu: float, seed_base: int) -> Dict[str, float]:
    sigma = max(float(np.std(H)), 1e-6)
    z0 = evidence_flat(H, sigma, n_mc=10000, seed=seed_base + 1, like_mode=like_mode, nu=nu)
    z1 = evidence_legacy(helper, theta, H, sigma, n_mc=13000, seed=seed_base + 2, like_mode=like_mode, nu=nu)
    z2 = evidence_reparam(helper, theta, H, sigma, q_center=q_center, q_width=q_width, n_mc=17000, seed=seed_base + 3, like_mode=like_mode, nu=nu)
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
    n_rep = 16
    nu = 4.0

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
            rows.append({"theta_deg": float(sep), "hxy": float(hxy)})

    if len(rows) < 85:
        raise RuntimeError("Base operational cohort too small.")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)

    full_gauss = run_family(helper, theta_all, H_all, q_center=q_center, q_width=q_width, like_mode="gauss", nu=nu, seed_base=10000)
    full_t = run_family(helper, theta_all, H_all, q_center=q_center, q_width=q_width, like_mode="student_t", nu=nu, seed_base=11000)

    rng = np.random.default_rng(1791)
    rep_rows = []
    for i in range(n_rep):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 80:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        r_g = run_family(helper, th, hh, q_center=q_center, q_width=q_width, like_mode="gauss", nu=nu, seed_base=12000 + 20 * i)
        r_t = run_family(helper, th, hh, q_center=q_center, q_width=q_width, like_mode="student_t", nu=nu, seed_base=13000 + 20 * i)
        rep_rows.append(
            {
                "rep": i,
                "n_pairs": int(len(idx)),
                "gauss_logB_reparam_vs_flat": r_g["logB_reparam_vs_flat"],
                "t_logB_reparam_vs_flat": r_t["logB_reparam_vs_flat"],
                "delta_t_minus_gauss": float(r_t["logB_reparam_vs_flat"] - r_g["logB_reparam_vs_flat"]),
            }
        )

    arr_t = np.array([r["t_logB_reparam_vs_flat"] for r in rep_rows], dtype=float)
    arr_d = np.array([r["delta_t_minus_gauss"] for r in rep_rows], dtype=float)

    summary = {
        "n_pairs_operational": len(rows),
        "q_width": q_width,
        "fraction": frac,
        "student_t_nu": nu,
        "replications": len(rep_rows),
        "full_gauss_logB_reparam_vs_flat": full_gauss["logB_reparam_vs_flat"],
        "full_t_logB_reparam_vs_flat": full_t["logB_reparam_vs_flat"],
        "full_delta_t_minus_gauss": float(full_t["logB_reparam_vs_flat"] - full_gauss["logB_reparam_vs_flat"]),
        "full_t_delta_reparam_minus_legacy": full_t["delta_logB_reparam_minus_legacy"],
        "prob_t_positive": float(np.mean(arr_t > 0.0)),
        "prob_t_better_than_gauss": float(np.mean(arr_d > 0.0)),
        "median_delta_t_minus_gauss": float(np.median(arr_d)),
        "std_delta_t_minus_gauss": float(np.std(arr_d)),
    }

    pass_full = summary["full_t_logB_reparam_vs_flat"] > 0.0
    pass_improvement = summary["full_delta_t_minus_gauss"] > 0.0
    pass_rep_positive = summary["prob_t_positive"] >= 0.95
    pass_rep_improve = summary["prob_t_better_than_gauss"] >= 0.75
    pass_disp = summary["std_delta_t_minus_gauss"] <= 0.15

    if pass_full and pass_improvement and pass_rep_positive and pass_rep_improve and pass_disp:
        verdict = "ROBUST_LIKELIHOOD_SUPPORTED"
    elif pass_full and pass_rep_positive:
        verdict = "ROBUST_LIKELIHOOD_PARTIAL"
    else:
        verdict = "ROBUST_LIKELIHOOD_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "pass_flags": {
            "full_positive": bool(pass_full),
            "full_improvement_over_gauss": bool(pass_improvement),
            "rep_positive": bool(pass_rep_positive),
            "rep_improvement": bool(pass_rep_improve),
            "dispersion_control": bool(pass_disp),
        },
        "verdict": verdict,
        "full_gaussian": full_gauss,
        "full_student_t": full_t,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1791: ROBUST LIKELIHOOD REPARAM",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Operational pairs: {len(rows)}",
        f"- Full gauss logB(reparam vs flat): {summary['full_gauss_logB_reparam_vs_flat']:.4f}",
        f"- Full Student-t logB(reparam vs flat): {summary['full_t_logB_reparam_vs_flat']:.4f}",
        f"- Full delta (t-gauss): {summary['full_delta_t_minus_gauss']:.4f}",
        f"- P(t>0): {summary['prob_t_positive']:.3f}",
        f"- P(t>gauss): {summary['prob_t_better_than_gauss']:.3f}",
        f"- Std delta (t-gauss): {summary['std_delta_t_minus_gauss']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_positive: {pass_full}",
        f"- full_improvement_over_gauss: {pass_improvement}",
        f"- rep_positive: {pass_rep_positive}",
        f"- rep_improvement: {pass_rep_improve}",
        f"- dispersion_control: {pass_disp}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1791] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1791] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
