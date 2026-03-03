#!/usr/bin/env python3
"""
QW-1811: Dynamic triad model (redesigned phase-1 representation).

Dynamic features from split-halves:
  f_drift      = h2 - h1
  f_persist    = sign(h1*h2) * sqrt(|h1*h2|)
  f_volatility = |h2 - h1|

Model:
  M2T(theta, f) = A * HD(theta)^q + b1*f_drift + b2*f_persist + b3*f_volatility + C
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
OUT_JSON = ROOT / "report_qw1811_dynamic_triad_model.json"
OUT_MD = ROOT / "RAPORT_QW1811_DYNAMIC_TRIAD_MODEL.md"


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


def center(v: np.ndarray) -> np.ndarray:
    return v - float(np.mean(v))


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


def evidence_m2t(helper, theta: np.ndarray, H: np.ndarray, F1: np.ndarray, F2: np.ndarray, F3: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B1 = rng.uniform(-1.0, 1.0, n_mc)
    B2 = rng.uniform(-1.0, 1.0, n_mc)
    B3 = rng.uniform(-1.0, 1.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array(
        [loglike(H, sigma, a * (hd0 ** qq) + b1 * F1 + b2 * F2 + b3 * F3 + c) for a, b1, b2, b3, c, qq in zip(A, B1, B2, B3, C, q)],
        dtype=float,
    )
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


def fit_best_m2t(helper, theta: np.ndarray, H: np.ndarray, F1: np.ndarray, F2: np.ndarray, F3: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B1 = rng.uniform(-1.0, 1.0, n_mc)
    B2 = rng.uniform(-1.0, 1.0, n_mc)
    B3 = rng.uniform(-1.0, 1.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    best_ll = -np.inf
    best = None
    for a, b1, b2, b3, c, qq in zip(A, B1, B2, B3, C, q):
        m = a * (hd0 ** qq) + b1 * F1 + b2 * F2 + b3 * F3 + c
        lk = loglike(H, sigma, m)
        if lk > best_ll:
            best_ll = lk
            best = {"A": float(a), "B1": float(b1), "B2": float(b2), "B3": float(b3), "C": float(c), "q": float(qq), "sigma": sigma}
    return best


def max_abs_feature_corr(resid: np.ndarray, features: List[np.ndarray]) -> float:
    vals = [abs(float(spearmanr(resid, f)[0])) for f in features]
    return float(max(vals))


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
            persist = float(np.sign(h1 * h2) * np.sqrt(abs(h1 * h2)))
            rows.append(
                {
                    "theta_deg": float(sep),
                    "hxy": float(hxy),
                    "f_drift": drift,
                    "f_persist": persist,
                    "f_vol": float(abs(drift)),
                }
            )

    if len(rows) < 85:
        raise RuntimeError("Operational cohort too small.")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    F1_all = center(np.array([r["f_drift"] for r in rows], dtype=float))
    F2_all = center(np.array([r["f_persist"] for r in rows], dtype=float))
    F3_all = center(np.array([r["f_vol"] for r in rows], dtype=float))
    hd0 = np.clip(helper.hellings_downs(theta_all), 1e-9, None)

    z0 = evidence_flat(H_all, n_mc=8500, seed=19101)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=19102)
    z2t = evidence_m2t(helper, theta_all, H_all, F1_all, F2_all, F3_all, q_center=q_center, q_width=q_width, n_mc=19000, seed=19103)
    full_m2 = float(z2 - z0)
    full_m2t = float(z2t - z0)
    full_delta = float(full_m2t - full_m2)

    fit2 = fit_best_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=65000, seed=19104)
    fit2t = fit_best_m2t(helper, theta_all, H_all, F1_all, F2_all, F3_all, q_center=q_center, q_width=q_width, n_mc=95000, seed=19105)
    pred2 = fit2["A"] * (hd0 ** fit2["q"]) + fit2["C"]
    pred2t = fit2t["A"] * (hd0 ** fit2t["q"]) + fit2t["B1"] * F1_all + fit2t["B2"] * F2_all + fit2t["B3"] * F3_all + fit2t["C"]
    resid2 = H_all - pred2
    resid2t = H_all - pred2t
    corr2 = max_abs_feature_corr(resid2, [F1_all, F2_all, F3_all])
    corr2t = max_abs_feature_corr(resid2t, [F1_all, F2_all, F3_all])
    feature_corr_improvement = float(corr2 - corr2t)

    rng = np.random.default_rng(19106)
    rep_rows = []
    for i in range(12):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 80:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        f1 = F1_all[idx]
        f2 = F2_all[idx]
        f3 = F3_all[idx]
        b0 = evidence_flat(hh, n_mc=3900, seed=19110 + 20 * i + 1)
        b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=19110 + 20 * i + 2)
        bt = evidence_m2t(helper, th, hh, f1, f2, f3, q_center=q_center, q_width=q_width, n_mc=8600, seed=19110 + 20 * i + 3)
        l2 = float(b2 - b0)
        lt = float(bt - b0)
        rep_rows.append({"rep": i, "logB_m2_vs_flat": l2, "logB_m2t_vs_flat": lt, "delta_m2t_vs_m2": float(lt - l2)})

    arr_d = np.array([r["delta_m2t_vs_m2"] for r in rep_rows], dtype=float)
    arr_t = np.array([r["logB_m2t_vs_flat"] for r in rep_rows], dtype=float)

    summary = {
        "n_pairs": len(rows),
        "fraction": frac,
        "q_width": q_width,
        "full_logB_m2_vs_flat": full_m2,
        "full_logB_m2t_vs_flat": full_m2t,
        "full_delta_m2t_vs_m2": full_delta,
        "replications": len(rep_rows),
        "prob_m2t_gt_m2": float(np.mean(arr_d > 0.0)),
        "prob_m2t_gt_flat": float(np.mean(arr_t > 0.0)),
        "median_delta_m2t_vs_m2": float(np.median(arr_d)),
        "std_delta_m2t_vs_m2": float(np.std(arr_d)),
        "max_abs_feature_corr_resid_m2": corr2,
        "max_abs_feature_corr_resid_m2t": corr2t,
        "feature_corr_improvement": feature_corr_improvement,
    }

    pass_full = summary["full_delta_m2t_vs_m2"] > 0.0
    pass_rep = summary["prob_m2t_gt_m2"] >= 0.80 and summary["prob_m2t_gt_flat"] >= 0.95
    pass_disp = summary["std_delta_m2t_vs_m2"] <= 0.30
    pass_dyn = summary["feature_corr_improvement"] > 0.05

    if pass_full and pass_rep and pass_disp and pass_dyn:
        verdict = "DYNAMIC_TRIAD_MODEL_SUPPORTED"
    elif pass_full and pass_rep and pass_dyn:
        verdict = "DYNAMIC_TRIAD_MODEL_PARTIAL"
    else:
        verdict = "DYNAMIC_TRIAD_MODEL_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "fit_params": {"m2": fit2, "m2t": fit2t},
        "pass_flags": {
            "full_gain": bool(pass_full),
            "replication_gain": bool(pass_rep),
            "dispersion_control": bool(pass_disp),
            "dynamic_feature_residual_reduction": bool(pass_dyn),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1811: DYNAMIC TRIAD MODEL",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Full logB M2/M2T: {full_m2:.4f} / {full_m2t:.4f}",
        f"- Full delta M2T-M2: {full_delta:.4f}",
        f"- P(M2T>M2): {summary['prob_m2t_gt_m2']:.3f}",
        f"- P(M2T>flat): {summary['prob_m2t_gt_flat']:.3f}",
        f"- Std delta M2T-M2: {summary['std_delta_m2t_vs_m2']:.3f}",
        f"- Feature-corr improvement: {feature_corr_improvement:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_gain: {pass_full}",
        f"- replication_gain: {pass_rep}",
        f"- dispersion_control: {pass_disp}",
        f"- dynamic_feature_residual_reduction: {pass_dyn}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1811] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1811] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
