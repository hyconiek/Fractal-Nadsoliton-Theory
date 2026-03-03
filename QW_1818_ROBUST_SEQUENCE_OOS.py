#!/usr/bin/env python3
"""
QW-1818: Robust OOS validation for sequence embedding branch.

Builds on QW-1817 by adding robust preprocessing:
- train-only winsorization of embedding features,
- stronger ridge regularization grid,
- same repeated stratified OOS protocol.

Goal:
- keep predictive gain,
- reduce cross-split dispersion in test delta log-likelihood.
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1818_robust_sequence_oos.json"
OUT_MD = ROOT / "RAPORT_QW1818_ROBUST_SEQUENCE_OOS.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def winsor_standardize_train_test(
    E_train: np.ndarray,
    E_test: np.ndarray,
    q_low: float = 0.05,
    q_high: float = 0.95,
) -> Tuple[np.ndarray, np.ndarray]:
    lo = np.quantile(E_train, q_low, axis=0)
    hi = np.quantile(E_train, q_high, axis=0)

    Et = np.clip(E_train, lo, hi)
    Ev = np.clip(E_test, lo, hi)

    mu = np.mean(Et, axis=0, keepdims=True)
    sd = np.std(Et, axis=0, keepdims=True)
    sd = np.where(sd <= 1e-12, 1.0, sd)
    return (Et - mu) / sd, (Ev - mu) / sd


def best_q_alpha_m2e_robust(
    m1817,
    hd_train: np.ndarray,
    y_train: np.ndarray,
    E_train: np.ndarray,
    q_grid: np.ndarray,
    alpha_grid: np.ndarray,
    rng: np.random.Generator,
) -> Tuple[float, float, np.ndarray, float]:
    idx = np.arange(len(y_train))
    rng.shuffle(idx)
    k_val = max(12, int(round(0.22 * len(idx))))
    k_val = min(k_val, len(idx) - 12)
    val_idx = np.sort(idx[:k_val])
    tr_idx = np.sort(idx[k_val:])

    if len(tr_idx) < 20 or len(val_idx) < 10:
        tr_idx = np.arange(len(y_train))
        val_idx = np.arange(len(y_train))

    hd_tr, hd_val = hd_train[tr_idx], hd_train[val_idx]
    y_tr, y_val = y_train[tr_idx], y_train[val_idx]
    E_tr_raw, E_val_raw = E_train[tr_idx], E_train[val_idx]
    E_tr, E_val = winsor_standardize_train_test(E_tr_raw, E_val_raw)

    best = None
    for q in q_grid:
        x_tr = hd_tr ** q
        x_val = hd_val ** q
        X_tr = np.column_stack([x_tr, E_tr, np.ones_like(x_tr)])
        X_val = np.column_stack([x_val, E_val, np.ones_like(x_val)])

        p = X_tr.shape[1]
        mask = np.zeros(p)
        mask[1 : 1 + E_tr.shape[1]] = 1.0

        for alpha in alpha_grid:
            beta = m1817.fit_ridge(X_tr, y_tr, alpha=alpha, penalize_mask=mask)
            pred_tr = X_tr @ beta
            pred_val = X_val @ beta
            sigma_tr = max(float(np.std(y_tr - pred_tr)), 1e-6)
            ll_val = m1817.gaussian_loglike(y_val, pred_val, sigma_tr)
            if best is None or ll_val > best[0]:
                best = (ll_val, float(q), float(alpha))

    _, q_best, alpha_best = best

    x_full = hd_train ** q_best
    E_full, _ = winsor_standardize_train_test(E_train, E_train)
    X_full = np.column_stack([x_full, E_full, np.ones_like(x_full)])

    p = X_full.shape[1]
    mask = np.zeros(p)
    mask[1 : 1 + E_full.shape[1]] = 1.0
    beta_full = m1817.fit_ridge(X_full, y_train, alpha=alpha_best, penalize_mask=mask)
    sigma_full = max(float(np.std(y_train - (X_full @ beta_full))), 1e-6)
    return q_best, alpha_best, beta_full, sigma_full


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
        raise RuntimeError(f"Robust OOS cohort too small: {len(rows)}")

    theta = np.array([r["theta_deg"] for r in rows], dtype=float)
    y = np.array([r["hxy"] for r in rows], dtype=float)
    hd = np.clip(helper.hellings_downs(theta), 1e-9, None)

    feature_names = ["f_mean", "f_std", "f_slope", "f_quad", "f_spread", "f_autoc1", "f_switch"]
    E_raw = np.column_stack([np.array([r[k] for r in rows], dtype=float) for k in feature_names])

    q_grid = np.linspace(max(0.9, q_center - 0.25), min(2.3, q_center + 0.25), 17)
    alpha_grid = np.array([0.30, 0.60, 1.20, 2.40, 4.80, 9.60, 14.40], dtype=float)

    rng = np.random.default_rng(19601)
    n_rep = 14
    rep_rows = []

    for i in range(n_rep):
        tr_idx, te_idx = m1817.stratified_split_indices(theta, train_frac=0.72, rng=rng, n_bins=8)
        if len(tr_idx) < 55 or len(te_idx) < 20:
            continue

        hd_tr, hd_te = hd[tr_idx], hd[te_idx]
        y_tr, y_te = y[tr_idx], y[te_idx]
        E_tr_raw, E_te_raw = E_raw[tr_idx], E_raw[te_idx]

        # M2 baseline
        q2, beta2, sigma2 = m1817.best_q_m2(hd_tr, y_tr, q_grid=q_grid)
        X2_te = np.column_stack([hd_te ** q2, np.ones_like(hd_te)])
        pred2_te = X2_te @ beta2
        ll2_te = m1817.gaussian_loglike(y_te, pred2_te, sigma2)
        rmse2 = m1817.rmse(y_te, pred2_te)

        # Robust M2E
        qE, alphaE, betaE, sigmaE = best_q_alpha_m2e_robust(
            m1817,
            hd_train=hd_tr,
            y_train=y_tr,
            E_train=E_tr_raw,
            q_grid=q_grid,
            alpha_grid=alpha_grid,
            rng=rng,
        )
        E_tr_w, E_te_w = winsor_standardize_train_test(E_tr_raw, E_te_raw)
        XE_te = np.column_stack([hd_te ** qE, E_te_w, np.ones_like(hd_te)])
        predE_te = XE_te @ betaE
        llE_te = m1817.gaussian_loglike(y_te, predE_te, sigmaE)
        rmseE = m1817.rmse(y_te, predE_te)

        # Flat
        c0 = float(np.mean(y_tr))
        sigma0 = max(float(np.std(y_tr - c0)), 1e-6)
        pred0_te = np.full_like(y_te, c0)
        ll0_te = m1817.gaussian_loglike(y_te, pred0_te, sigma0)

        rep_rows.append(
            {
                "rep": i,
                "n_train": int(len(tr_idx)),
                "n_test": int(len(te_idx)),
                "q_m2": float(q2),
                "q_m2e": float(qE),
                "alpha_m2e": float(alphaE),
                "ll_test_flat": float(ll0_te),
                "ll_test_m2": float(ll2_te),
                "ll_test_m2e": float(llE_te),
                "delta_ll_m2_vs_flat": float(ll2_te - ll0_te),
                "delta_ll_m2e_vs_flat": float(llE_te - ll0_te),
                "delta_ll_m2e_vs_m2": float(llE_te - ll2_te),
                "rmse_m2": float(rmse2),
                "rmse_m2e": float(rmseE),
                "rmse_gain_m2_minus_m2e": float(rmse2 - rmseE),
            }
        )

    if len(rep_rows) < 10:
        raise RuntimeError("Too few robust OOS replications.")

    arr_delta = np.array([r["delta_ll_m2e_vs_m2"] for r in rep_rows], dtype=float)
    arr_gain = np.array([r["rmse_gain_m2_minus_m2e"] for r in rep_rows], dtype=float)
    arr_flat = np.array([r["delta_ll_m2e_vs_flat"] for r in rep_rows], dtype=float)

    std_prev = float(d1817["summary"]["std_delta_ll_m2e_vs_m2"])
    std_now = float(np.std(arr_delta))

    summary = {
        "n_pairs_oos_cohort": len(rows),
        "mean_windows_per_pair": float(np.mean([r["n_windows"] for r in rows])),
        "n_features": int(E_raw.shape[1]),
        "n_replications": len(rep_rows),
        "train_fraction": 0.72,
        "mean_delta_ll_m2e_vs_m2": float(np.mean(arr_delta)),
        "median_delta_ll_m2e_vs_m2": float(np.median(arr_delta)),
        "std_delta_ll_m2e_vs_m2": std_now,
        "prob_delta_ll_m2e_gt_m2": float(np.mean(arr_delta > 0.0)),
        "prob_m2e_gt_flat": float(np.mean(arr_flat > 0.0)),
        "mean_rmse_gain_m2_minus_m2e": float(np.mean(arr_gain)),
        "prob_rmse_gain_positive": float(np.mean(arr_gain > 0.0)),
        "std_prev_qw1817": std_prev,
        "std_reduction_vs_qw1817": float(std_prev - std_now),
        "std_reduction_ratio_vs_qw1817": float((std_prev - std_now) / max(std_prev, 1e-12)),
    }

    pass_ll = summary["mean_delta_ll_m2e_vs_m2"] > 0.0 and summary["prob_delta_ll_m2e_gt_m2"] >= 0.75
    pass_rmse = summary["mean_rmse_gain_m2_minus_m2e"] > 0.0 and summary["prob_rmse_gain_positive"] >= 0.75
    pass_disp = summary["std_reduction_ratio_vs_qw1817"] >= 0.25

    if pass_ll and pass_rmse and pass_disp:
        verdict = "ROBUST_SEQUENCE_OOS_SUPPORTED"
    elif pass_ll and pass_rmse:
        verdict = "ROBUST_SEQUENCE_OOS_PARTIAL"
    else:
        verdict = "ROBUST_SEQUENCE_OOS_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "feature_names": feature_names,
        "summary": summary,
        "pass_flags": {
            "test_ll_gain": bool(pass_ll),
            "test_rmse_gain": bool(pass_rmse),
            "dispersion_reduction_vs_qw1817": bool(pass_disp),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1818: ROBUST SEQUENCE OOS",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- OOS cohort size: {len(rows)}",
        f"- Mean windows per pair: {summary['mean_windows_per_pair']:.2f}",
        f"- Features: {', '.join(feature_names)}",
        f"- Replications: {summary['n_replications']}",
        f"- Mean test delta LL (M2E-M2): {summary['mean_delta_ll_m2e_vs_m2']:.4f}",
        f"- P(test LL M2E>M2): {summary['prob_delta_ll_m2e_gt_m2']:.3f}",
        f"- Mean RMSE gain (M2-M2E): {summary['mean_rmse_gain_m2_minus_m2e']:.6f}",
        f"- P(RMSE gain>0): {summary['prob_rmse_gain_positive']:.3f}",
        f"- Std delta LL (QW-1817 -> QW-1818): {summary['std_prev_qw1817']:.3f} -> {summary['std_delta_ll_m2e_vs_m2']:.3f}",
        f"- Std reduction ratio: {summary['std_reduction_ratio_vs_qw1817']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- test_ll_gain: {pass_ll}",
        f"- test_rmse_gain: {pass_rmse}",
        f"- dispersion_reduction_vs_qw1817: {pass_disp}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1818] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1818] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
