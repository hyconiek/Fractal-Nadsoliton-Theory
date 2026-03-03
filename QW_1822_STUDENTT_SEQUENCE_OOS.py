#!/usr/bin/env python3
"""
QW-1822: Student-t OOS scoring for sequence embedding branch.

Implements recommendation from QW-1821:
- keep mean models from QW-1817,
- replace Gaussian test score with Student-t predictive log-score,
- evaluate whether heavy tails reduce split dispersion.
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
from scipy.special import gammaln


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1822_studentt_sequence_oos.json"
OUT_MD = ROOT / "RAPORT_QW1822_STUDENTT_SEQUENCE_OOS.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def mad_scale(resid: np.ndarray) -> float:
    med = np.median(resid)
    mad = np.median(np.abs(resid - med))
    s = 1.4826 * mad
    if not np.isfinite(s) or s <= 1e-8:
        s = float(np.std(resid))
    return max(float(s), 1e-6)


def student_t_loglike(y: np.ndarray, mu: np.ndarray, sigma: float, nu: float) -> float:
    z = (y - mu) / sigma
    n = len(y)
    c = gammaln((nu + 1.0) / 2.0) - gammaln(nu / 2.0) - 0.5 * np.log(np.pi * nu) - np.log(sigma)
    return float(np.sum(c - ((nu + 1.0) / 2.0) * np.log1p((z * z) / nu)))


def best_nu_on_train(y: np.ndarray, mu: np.ndarray, sigma: float, nu_grid: np.ndarray) -> float:
    best = None
    for nu in nu_grid:
        ll = student_t_loglike(y, mu, sigma, float(nu))
        if best is None or ll > best[0]:
            best = (ll, float(nu))
    return best[1]


def main() -> None:
    helper = load_module(ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py", "qw1787_helper")
    m1817 = load_module(ROOT / "QW_1817_SEQUENCE_OOS_VALIDATION.py", "qw1817_oos")

    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))
    d1817 = json.loads((ROOT / "report_qw1817_sequence_oos_validation.json").read_text(encoding="utf-8"))

    q_center = float(d1773["projection"]["p"])
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
        if len(x) < n_match_min or stab > stability_max:
            continue

        t, v = m1817.windowed_multiscale_traj(helper, x, y)
        feats = m1817.trajectory_features(t, v)
        if feats is None:
            continue

        rows.append(
            {
                "theta_deg": float(sep),
                "hxy": float(hxy),
                "f_mean": feats["f_mean"],
                "f_std": feats["f_std"],
                "f_slope": feats["f_slope"],
                "f_quad": feats["f_quad"],
                "f_spread": feats["f_spread"],
                "f_autoc1": feats["f_autoc1"],
                "f_switch": feats["f_switch"],
                "n_windows": feats["n_windows"],
            }
        )

    if len(rows) < 80:
        raise RuntimeError(f"Student-t OOS cohort too small: {len(rows)}")

    theta = np.array([r["theta_deg"] for r in rows], dtype=float)
    y = np.array([r["hxy"] for r in rows], dtype=float)
    hd = np.clip(helper.hellings_downs(theta), 1e-9, None)

    feature_names = ["f_mean", "f_std", "f_slope", "f_quad", "f_spread", "f_autoc1", "f_switch"]
    E_raw = np.column_stack([np.array([r[k] for r in rows], dtype=float) for k in feature_names])

    q_grid = np.linspace(max(0.9, q_center - 0.25), min(2.3, q_center + 0.25), 17)
    alpha_grid = np.array([0.08, 0.15, 0.30, 0.60, 1.20, 2.40], dtype=float)
    nu_grid = np.array([3.0, 4.0, 5.0, 8.0, 12.0, 20.0, 40.0], dtype=float)

    rng = np.random.default_rng(19801)
    n_rep = 14
    rep_rows = []

    for i in range(n_rep):
        tr_idx, te_idx = m1817.stratified_split_indices(theta, train_frac=0.72, rng=rng, n_bins=8)
        if len(tr_idx) < 55 or len(te_idx) < 20:
            continue

        hd_tr, hd_te = hd[tr_idx], hd[te_idx]
        y_tr, y_te = y[tr_idx], y[te_idx]
        E_tr_raw, E_te_raw = E_raw[tr_idx], E_raw[te_idx]
        E_tr, E_te = m1817.standardize_train_test(E_tr_raw, E_te_raw)

        # M2 mean model
        q2, beta2, _ = m1817.best_q_m2(hd_tr, y_tr, q_grid=q_grid)
        X2_tr = np.column_stack([hd_tr ** q2, np.ones_like(hd_tr)])
        X2_te = np.column_stack([hd_te ** q2, np.ones_like(hd_te)])
        pred2_tr = X2_tr @ beta2
        pred2_te = X2_te @ beta2
        s2 = mad_scale(y_tr - pred2_tr)
        nu2 = best_nu_on_train(y_tr, pred2_tr, s2, nu_grid=nu_grid)
        ll2_te = student_t_loglike(y_te, pred2_te, s2, nu2)

        # M2E mean model
        qE, alphaE, betaE, sigmaE = m1817.best_q_alpha_m2e(hd_tr, y_tr, E_tr, q_grid=q_grid, alpha_grid=alpha_grid, rng=rng)
        XE_tr = np.column_stack([hd_tr ** qE, E_tr, np.ones_like(hd_tr)])
        XE_te = np.column_stack([hd_te ** qE, E_te, np.ones_like(hd_te)])
        predE_tr = XE_tr @ betaE
        predE_te = XE_te @ betaE
        sE = mad_scale(y_tr - predE_tr)
        nuE = best_nu_on_train(y_tr, predE_tr, sE, nu_grid=nu_grid)
        llE_te = student_t_loglike(y_te, predE_te, sE, nuE)

        # Flat baseline
        c0 = float(np.mean(y_tr))
        pred0_tr = np.full_like(y_tr, c0)
        pred0_te = np.full_like(y_te, c0)
        s0 = mad_scale(y_tr - pred0_tr)
        nu0 = best_nu_on_train(y_tr, pred0_tr, s0, nu_grid=nu_grid)
        ll0_te = student_t_loglike(y_te, pred0_te, s0, nu0)

        rep_rows.append(
            {
                "rep": i,
                "n_train": int(len(tr_idx)),
                "n_test": int(len(te_idx)),
                "q_m2": float(q2),
                "q_m2e": float(qE),
                "alpha_m2e": float(alphaE),
                "nu_flat": float(nu0),
                "nu_m2": float(nu2),
                "nu_m2e": float(nuE),
                "ll_t_test_flat": float(ll0_te),
                "ll_t_test_m2": float(ll2_te),
                "ll_t_test_m2e": float(llE_te),
                "delta_ll_t_m2e_vs_m2": float(llE_te - ll2_te),
                "delta_ll_t_m2e_vs_flat": float(llE_te - ll0_te),
                "delta_ll_t_m2_vs_flat": float(ll2_te - ll0_te),
                "sigma_m2e_train": float(sE),
                "sigma_m2e_gauss_train": float(sigmaE),
            }
        )

    if len(rep_rows) < 10:
        raise RuntimeError("Too few Student-t OOS replications.")

    arr_delta = np.array([r["delta_ll_t_m2e_vs_m2"] for r in rep_rows], dtype=float)
    arr_flat = np.array([r["delta_ll_t_m2e_vs_flat"] for r in rep_rows], dtype=float)

    std_prev = float(d1817["summary"]["std_delta_ll_m2e_vs_m2"])
    std_now = float(np.std(arr_delta))

    summary = {
        "n_pairs_oos_cohort": len(rows),
        "mean_windows_per_pair": float(np.mean([r["n_windows"] for r in rows])),
        "n_features": int(E_raw.shape[1]),
        "n_replications": len(rep_rows),
        "train_fraction": 0.72,
        "mean_delta_ll_t_m2e_vs_m2": float(np.mean(arr_delta)),
        "median_delta_ll_t_m2e_vs_m2": float(np.median(arr_delta)),
        "std_delta_ll_t_m2e_vs_m2": std_now,
        "prob_delta_ll_t_positive": float(np.mean(arr_delta > 0.0)),
        "prob_m2e_t_gt_flat": float(np.mean(arr_flat > 0.0)),
        "mean_nu_m2e": float(np.mean([r["nu_m2e"] for r in rep_rows])),
        "std_prev_qw1817": std_prev,
        "std_reduction_vs_qw1817": float(std_prev - std_now),
        "std_reduction_ratio_vs_qw1817": float((std_prev - std_now) / max(std_prev, 1e-12)),
    }

    pass_gain = summary["mean_delta_ll_t_m2e_vs_m2"] > 0.0 and summary["prob_delta_ll_t_positive"] >= 0.80
    pass_flat = summary["prob_m2e_t_gt_flat"] >= 0.90
    pass_disp = summary["std_reduction_ratio_vs_qw1817"] >= 0.20

    if pass_gain and pass_flat and pass_disp:
        verdict = "STUDENTT_SEQUENCE_OOS_SUPPORTED"
    elif pass_gain and pass_flat:
        verdict = "STUDENTT_SEQUENCE_OOS_PARTIAL"
    else:
        verdict = "STUDENTT_SEQUENCE_OOS_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "feature_names": feature_names,
        "summary": summary,
        "pass_flags": {
            "studentt_gain_vs_m2": bool(pass_gain),
            "studentt_gain_vs_flat": bool(pass_flat),
            "dispersion_reduction_vs_qw1817": bool(pass_disp),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1822: STUDENT-T SEQUENCE OOS",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- OOS cohort size: {len(rows)}",
        f"- Mean windows per pair: {summary['mean_windows_per_pair']:.2f}",
        f"- Features: {', '.join(feature_names)}",
        f"- Replications: {summary['n_replications']}",
        f"- Mean delta t-LL (M2E-M2): {summary['mean_delta_ll_t_m2e_vs_m2']:.4f}",
        f"- P(delta t-LL > 0): {summary['prob_delta_ll_t_positive']:.3f}",
        f"- P(M2E_t > flat_t): {summary['prob_m2e_t_gt_flat']:.3f}",
        f"- Mean nu(M2E): {summary['mean_nu_m2e']:.3f}",
        f"- Std delta LL (QW-1817 -> Student-t): {summary['std_prev_qw1817']:.3f} -> {summary['std_delta_ll_t_m2e_vs_m2']:.3f}",
        f"- Std reduction ratio: {summary['std_reduction_ratio_vs_qw1817']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- studentt_gain_vs_m2: {pass_gain}",
        f"- studentt_gain_vs_flat: {pass_flat}",
        f"- dispersion_reduction_vs_qw1817: {pass_disp}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1822] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1822] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
