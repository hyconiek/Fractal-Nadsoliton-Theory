#!/usr/bin/env python3
"""
QW-1806: Quality-regime diagnostic model.

Hypothesis:
- high split-to-split instability (QW-1804) may be driven by hidden composition shifts
  in pair quality (n_match, split-half stability), not only angular physics.

Diagnostic model:
    M2Q(theta, Q) = A * HD(theta)^q + B * Q + C
where Q is centered quality proxy.
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
OUT_JSON = ROOT / "report_qw1806_quality_regime_diagnostic.json"
OUT_MD = ROOT / "RAPORT_QW1806_QUALITY_REGIME_DIAGNOSTIC.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
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
    q = 0.5 * (minmax(inv_sqrt_n) + minmax(stability))
    return q - float(np.mean(q))


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


def evidence_m2q(helper, theta: np.ndarray, H: np.ndarray, Q: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.0, 1.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array([loglike(H, sigma, a * (hd0 ** qq) + b * Q + c) for a, b, c, qq in zip(A, B, C, q)], dtype=float)
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


def fit_best_m2q(helper, theta: np.ndarray, H: np.ndarray, Q: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.0, 1.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    best_ll = -np.inf
    best = None
    for a, b, c, qq in zip(A, B, C, q):
        m = a * (hd0 ** qq) + b * Q + c
        lk = loglike(H, sigma, m)
        if lk > best_ll:
            best_ll = lk
            best = {"A": float(a), "B": float(b), "C": float(c), "q": float(qq), "sigma": sigma}
    return best


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

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    nm_all = np.array([r["n_match"] for r in rows], dtype=float)
    st_all = np.array([r["stability"] for r in rows], dtype=float)
    Q_all = quality_proxy(nm_all, st_all)
    hd0 = np.clip(helper.hellings_downs(theta_all), 1e-9, None)

    z0 = evidence_flat(H_all, n_mc=8500, seed=18601)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=18602)
    z2q = evidence_m2q(helper, theta_all, H_all, Q_all, q_center=q_center, q_width=q_width, n_mc=17000, seed=18603)
    full_m2 = float(z2 - z0)
    full_m2q = float(z2q - z0)
    full_delta = float(full_m2q - full_m2)

    fit2 = fit_best_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=65000, seed=18604)
    fit2q = fit_best_m2q(helper, theta_all, H_all, Q_all, q_center=q_center, q_width=q_width, n_mc=85000, seed=18605)
    pred2 = fit2["A"] * (hd0 ** fit2["q"]) + fit2["C"]
    pred2q = fit2q["A"] * (hd0 ** fit2q["q"]) + fit2q["B"] * Q_all + fit2q["C"]
    resid2 = H_all - pred2
    resid2q = H_all - pred2q
    rho_q_resid2 = float(spearmanr(Q_all, resid2)[0])
    rho_q_resid2q = float(spearmanr(Q_all, resid2q)[0])
    quality_improvement = float(abs(rho_q_resid2) - abs(rho_q_resid2q))

    rng = np.random.default_rng(18606)
    rep_rows = []
    for i in range(12):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 80:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        qq = Q_all[idx]
        b0 = evidence_flat(hh, n_mc=3900, seed=18610 + 20 * i + 1)
        b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=18610 + 20 * i + 2)
        b2q = evidence_m2q(helper, th, hh, qq, q_center=q_center, q_width=q_width, n_mc=8200, seed=18610 + 20 * i + 3)
        l2 = float(b2 - b0)
        l2q = float(b2q - b0)
        rep_rows.append({"rep": i, "logB_m2_vs_flat": l2, "logB_m2q_vs_flat": l2q, "delta_m2q_vs_m2": float(l2q - l2)})

    arr_q = np.array([r["logB_m2q_vs_flat"] for r in rep_rows], dtype=float)
    arr_d = np.array([r["delta_m2q_vs_m2"] for r in rep_rows], dtype=float)

    summary = {
        "n_pairs": len(rows),
        "fraction": frac,
        "q_width": q_width,
        "full_logB_m2_vs_flat": full_m2,
        "full_logB_m2q_vs_flat": full_m2q,
        "full_delta_m2q_vs_m2": full_delta,
        "replications": len(rep_rows),
        "prob_m2q_gt_m2": float(np.mean(arr_d > 0.0)),
        "prob_m2q_gt_flat": float(np.mean(arr_q > 0.0)),
        "median_delta_m2q_vs_m2": float(np.median(arr_d)),
        "std_delta_m2q_vs_m2": float(np.std(arr_d)),
        "rho_Q_resid_m2": rho_q_resid2,
        "rho_Q_resid_m2q": rho_q_resid2q,
        "quality_residual_improvement_abs": quality_improvement,
    }

    pass_full = summary["full_delta_m2q_vs_m2"] > 0.0
    pass_rep = summary["prob_m2q_gt_m2"] >= 0.80 and summary["prob_m2q_gt_flat"] >= 0.95
    pass_disp = summary["std_delta_m2q_vs_m2"] <= 0.30
    pass_quality = summary["quality_residual_improvement_abs"] > 0.05

    if pass_full and pass_rep and pass_disp and pass_quality:
        verdict = "QUALITY_REGIME_SIGNAL_SUPPORTED"
    elif pass_full and pass_rep and pass_quality:
        verdict = "QUALITY_REGIME_SIGNAL_PARTIAL"
    else:
        verdict = "QUALITY_REGIME_SIGNAL_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "fit_params": {"m2": fit2, "m2q": fit2q},
        "pass_flags": {
            "full_gain": bool(pass_full),
            "replication_gain": bool(pass_rep),
            "dispersion_control": bool(pass_disp),
            "quality_residual_reduction": bool(pass_quality),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1806: QUALITY REGIME DIAGNOSTIC",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Full logB M2/M2Q: {full_m2:.4f} / {full_m2q:.4f}",
        f"- Full delta M2Q-M2: {full_delta:.4f}",
        f"- P(M2Q>M2): {summary['prob_m2q_gt_m2']:.3f}",
        f"- P(M2Q>flat): {summary['prob_m2q_gt_flat']:.3f}",
        f"- Std delta M2Q-M2: {summary['std_delta_m2q_vs_m2']:.3f}",
        f"- Quality residual improvement |rho|: {quality_improvement:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_gain: {pass_full}",
        f"- replication_gain: {pass_rep}",
        f"- dispersion_control: {pass_disp}",
        f"- quality_residual_reduction: {pass_quality}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1806] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1806] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
