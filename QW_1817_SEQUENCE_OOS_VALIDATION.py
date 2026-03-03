#!/usr/bin/env python3
"""
QW-1817: Out-of-sample validation for sequence embedding branch.

Purpose:
- verify whether the strong in-sample Bayes evidence from QW-1815 translates
  into predictive out-of-sample gain,
- reduce risk of flexible-model overfitting.

Compared models:
- M2   : y = A * HD(theta)^q + C
- M2E  : y = A * HD(theta)^q + <B, E_std> + C  (ridge-regularized on B)
- Flat : y = C

Protocol:
- repeated stratified train/test splits over angular bins,
- inner validation on train for selecting (q, alpha) in M2E,
- test-set log-likelihood and RMSE comparison.
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
OUT_JSON = ROOT / "report_qw1817_sequence_oos_validation.json"
OUT_MD = ROOT / "RAPORT_QW1817_SEQUENCE_OOS_VALIDATION.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def windowed_multiscale_traj(helper, x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n = min(len(x), len(y))
    if n < 150:
        return np.array([], dtype=float), np.array([], dtype=float)

    scales = [0.45, 0.55, 0.65]
    centers = []
    vals = []
    for frac in scales:
        w = max(70, int(round(frac * n)))
        step = max(10, int(round(0.18 * w)))
        for start in range(0, n - w + 1, step):
            xx = x[start : start + w]
            yy = y[start : start + w]
            h = helper.cross_dfa(xx, yy, min_scale=10)
            if not np.isfinite(h):
                continue
            c = (start + 0.5 * w) / n
            centers.append(float(c))
            vals.append(float(h))
    if len(vals) < 8:
        return np.array([], dtype=float), np.array([], dtype=float)
    return np.array(centers, dtype=float), np.array(vals, dtype=float)


def lag1_autocorr(v: np.ndarray) -> float:
    if len(v) < 3:
        return 0.0
    a = v[:-1]
    b = v[1:]
    sa = float(np.std(a))
    sb = float(np.std(b))
    if sa <= 1e-12 or sb <= 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def switch_rate(v: np.ndarray) -> float:
    if len(v) < 2:
        return 0.0
    s = np.sign(v)
    return float(np.sum(s[1:] != s[:-1]) / max(len(v) - 1, 1))


def trajectory_features(centers: np.ndarray, vals: np.ndarray) -> Dict[str, float] | None:
    if len(vals) < 8:
        return None
    order = np.argsort(centers)
    t = centers[order]
    v = vals[order]

    p2 = np.polyfit(t, v, 2)
    quad = float(p2[0])
    slope = float(p2[1])
    mean = float(np.mean(v))
    std = float(np.std(v))
    p10, p90 = np.quantile(v, [0.10, 0.90])
    spread = float(p90 - p10)
    autoc = lag1_autocorr(v)
    sw = switch_rate(v - mean)

    return {
        "f_mean": mean,
        "f_std": std,
        "f_slope": slope,
        "f_quad": quad,
        "f_spread": spread,
        "f_autoc1": autoc,
        "f_switch": sw,
        "n_windows": int(len(v)),
    }


def stratified_split_indices(theta: np.ndarray, train_frac: float, rng: np.random.Generator, n_bins: int = 8) -> Tuple[np.ndarray, np.ndarray]:
    bins = np.linspace(0.0, 180.0, n_bins + 1)
    idx_all = np.arange(len(theta))
    tr = []
    te = []
    for i in range(n_bins):
        m = (theta >= bins[i]) & (theta < bins[i + 1] if i < n_bins - 1 else theta <= bins[i + 1])
        idx = idx_all[m]
        if len(idx) == 0:
            continue
        k = max(1, int(round(train_frac * len(idx))))
        k = min(k, len(idx) - 1) if len(idx) > 1 else 1
        idx_perm = rng.permutation(idx)
        tr.append(idx_perm[:k])
        te.append(idx_perm[k:])
    train = np.sort(np.concatenate([x for x in tr if len(x) > 0])) if tr else np.array([], dtype=int)
    test = np.sort(np.concatenate([x for x in te if len(x) > 0])) if te else np.array([], dtype=int)
    return train, test


def standardize_train_test(E_train: np.ndarray, E_test: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mu = np.mean(E_train, axis=0, keepdims=True)
    sd = np.std(E_train, axis=0, keepdims=True)
    sd = np.where(sd <= 1e-12, 1.0, sd)
    return (E_train - mu) / sd, (E_test - mu) / sd


def fit_ols(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


def fit_ridge(X: np.ndarray, y: np.ndarray, alpha: float, penalize_mask: np.ndarray) -> np.ndarray:
    xtx = X.T @ X
    reg = np.diag(penalize_mask.astype(float)) * alpha
    beta = np.linalg.solve(xtx + reg + 1e-9 * np.eye(xtx.shape[0]), X.T @ y)
    return beta


def rmse(y: np.ndarray, yhat: np.ndarray) -> float:
    return float(np.sqrt(np.mean((y - yhat) ** 2)))


def gaussian_loglike(y: np.ndarray, yhat: np.ndarray, sigma: float) -> float:
    z = (y - yhat) / sigma
    return float(-0.5 * np.sum(z * z) - len(y) * np.log(max(sigma, 1e-12)))


def best_q_m2(hd_train: np.ndarray, y_train: np.ndarray, q_grid: np.ndarray) -> Tuple[float, np.ndarray, float]:
    best = None
    for q in q_grid:
        x = hd_train ** q
        X = np.column_stack([x, np.ones_like(x)])
        beta = fit_ols(X, y_train)
        pred = X @ beta
        sigma = max(float(np.std(y_train - pred)), 1e-6)
        ll = gaussian_loglike(y_train, pred, sigma)
        if best is None or ll > best[0]:
            best = (ll, float(q), beta, sigma)
    _, q_best, beta_best, sigma_best = best
    return q_best, beta_best, sigma_best


def best_q_alpha_m2e(
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
    E_tr, E_val = E_train[tr_idx], E_train[val_idx]

    best = None
    for q in q_grid:
        x_tr = hd_tr ** q
        x_val = hd_val ** q
        X_tr = np.column_stack([x_tr, E_tr, np.ones_like(x_tr)])
        X_val = np.column_stack([x_val, E_val, np.ones_like(x_val)])

        # penalize only embedding coefficients
        p = X_tr.shape[1]
        mask = np.zeros(p)
        mask[1 : 1 + E_tr.shape[1]] = 1.0

        for alpha in alpha_grid:
            beta = fit_ridge(X_tr, y_tr, alpha=alpha, penalize_mask=mask)
            pred_val = X_val @ beta
            sigma_val = max(float(np.std(y_tr - (X_tr @ beta))), 1e-6)
            ll_val = gaussian_loglike(y_val, pred_val, sigma_val)
            if best is None or ll_val > best[0]:
                best = (ll_val, float(q), float(alpha))

    _, q_best, alpha_best = best

    x_full = hd_train ** q_best
    X_full = np.column_stack([x_full, E_train, np.ones_like(x_full)])
    p = X_full.shape[1]
    mask = np.zeros(p)
    mask[1 : 1 + E_train.shape[1]] = 1.0
    beta_full = fit_ridge(X_full, y_train, alpha=alpha_best, penalize_mask=mask)
    sigma_full = max(float(np.std(y_train - (X_full @ beta_full))), 1e-6)
    return q_best, alpha_best, beta_full, sigma_full


def main() -> None:
    helper = load_helper()
    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))

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

        t, v = windowed_multiscale_traj(helper, x, y)
        feats = trajectory_features(t, v)
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
        raise RuntimeError(f"OOS cohort too small: {len(rows)}")

    theta = np.array([r["theta_deg"] for r in rows], dtype=float)
    y = np.array([r["hxy"] for r in rows], dtype=float)
    hd = np.clip(helper.hellings_downs(theta), 1e-9, None)

    feature_names = ["f_mean", "f_std", "f_slope", "f_quad", "f_spread", "f_autoc1", "f_switch"]
    E_raw = np.column_stack([np.array([r[k] for r in rows], dtype=float) for k in feature_names])

    q_grid = np.linspace(max(0.9, q_center - 0.25), min(2.3, q_center + 0.25), 17)
    alpha_grid = np.array([0.08, 0.15, 0.30, 0.60, 1.20, 2.40], dtype=float)

    rng = np.random.default_rng(19501)
    n_rep = 14
    rep_rows = []

    for i in range(n_rep):
        tr_idx, te_idx = stratified_split_indices(theta, train_frac=0.72, rng=rng, n_bins=8)
        if len(tr_idx) < 55 or len(te_idx) < 20:
            continue

        hd_tr, hd_te = hd[tr_idx], hd[te_idx]
        y_tr, y_te = y[tr_idx], y[te_idx]
        E_tr_raw, E_te_raw = E_raw[tr_idx], E_raw[te_idx]
        E_tr, E_te = standardize_train_test(E_tr_raw, E_te_raw)

        # M2
        q2, beta2, sigma2 = best_q_m2(hd_tr, y_tr, q_grid=q_grid)
        X2_te = np.column_stack([hd_te ** q2, np.ones_like(hd_te)])
        pred2_te = X2_te @ beta2
        ll2_te = gaussian_loglike(y_te, pred2_te, sigma2)
        rmse2 = rmse(y_te, pred2_te)

        # M2E
        qE, alphaE, betaE, sigmaE = best_q_alpha_m2e(hd_tr, y_tr, E_tr, q_grid=q_grid, alpha_grid=alpha_grid, rng=rng)
        XE_te = np.column_stack([hd_te ** qE, E_te, np.ones_like(hd_te)])
        predE_te = XE_te @ betaE
        llE_te = gaussian_loglike(y_te, predE_te, sigmaE)
        rmseE = rmse(y_te, predE_te)

        # Flat
        c0 = float(np.mean(y_tr))
        sigma0 = max(float(np.std(y_tr - c0)), 1e-6)
        pred0_te = np.full_like(y_te, c0)
        ll0_te = gaussian_loglike(y_te, pred0_te, sigma0)

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
        raise RuntimeError("Too few valid OOS replications.")

    arr_delta = np.array([r["delta_ll_m2e_vs_m2"] for r in rep_rows], dtype=float)
    arr_gain = np.array([r["rmse_gain_m2_minus_m2e"] for r in rep_rows], dtype=float)
    arr_flat = np.array([r["delta_ll_m2e_vs_flat"] for r in rep_rows], dtype=float)

    summary = {
        "n_pairs_oos_cohort": len(rows),
        "mean_windows_per_pair": float(np.mean([r["n_windows"] for r in rows])),
        "n_features": int(E_raw.shape[1]),
        "n_replications": len(rep_rows),
        "train_fraction": 0.72,
        "mean_delta_ll_m2e_vs_m2": float(np.mean(arr_delta)),
        "median_delta_ll_m2e_vs_m2": float(np.median(arr_delta)),
        "std_delta_ll_m2e_vs_m2": float(np.std(arr_delta)),
        "prob_delta_ll_m2e_gt_m2": float(np.mean(arr_delta > 0.0)),
        "prob_m2e_gt_flat": float(np.mean(arr_flat > 0.0)),
        "mean_rmse_gain_m2_minus_m2e": float(np.mean(arr_gain)),
        "prob_rmse_gain_positive": float(np.mean(arr_gain > 0.0)),
    }

    pass_ll = summary["mean_delta_ll_m2e_vs_m2"] > 0.0 and summary["prob_delta_ll_m2e_gt_m2"] >= 0.75
    pass_rmse = summary["mean_rmse_gain_m2_minus_m2e"] > 0.0 and summary["prob_rmse_gain_positive"] >= 0.75
    pass_disp = summary["std_delta_ll_m2e_vs_m2"] <= 1.10

    if pass_ll and pass_rmse and pass_disp:
        verdict = "SEQUENCE_OOS_VALIDATION_SUPPORTED"
    elif pass_ll and pass_rmse:
        verdict = "SEQUENCE_OOS_VALIDATION_PARTIAL"
    else:
        verdict = "SEQUENCE_OOS_VALIDATION_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "feature_names": feature_names,
        "summary": summary,
        "pass_flags": {
            "test_ll_gain": bool(pass_ll),
            "test_rmse_gain": bool(pass_rmse),
            "dispersion_control": bool(pass_disp),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1817: SEQUENCE OOS VALIDATION",
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
        f"- Std test delta LL: {summary['std_delta_ll_m2e_vs_m2']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- test_ll_gain: {pass_ll}",
        f"- test_rmse_gain: {pass_rmse}",
        f"- dispersion_control: {pass_disp}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1817] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1817] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
