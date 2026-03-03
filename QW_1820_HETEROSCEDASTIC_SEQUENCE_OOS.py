#!/usr/bin/env python3
"""
QW-1820: Heteroscedastic OOS test for sequence embedding branch.

Redesign after QW-1819 inconsistency list:
- keep mean model M2E from QW-1817,
- replace constant-variance likelihood with feature-dependent variance model,
- test whether this improves test LL stability and reduces dispersion.
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
OUT_JSON = ROOT / "report_qw1820_heteroscedastic_sequence_oos.json"
OUT_MD = ROOT / "RAPORT_QW1820_HETEROSCEDASTIC_SEQUENCE_OOS.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fit_log_variance_model(resid: np.ndarray, E_std: np.ndarray) -> np.ndarray:
    """Fit log(res^2) ~ [1, |f_std|, |f_spread|, |f_switch|]."""
    eps = 1e-9
    z = np.log(resid * resid + eps)
    # feature indices from QW-1817 layout
    # [f_mean, f_std, f_slope, f_quad, f_spread, f_autoc1, f_switch]
    x_std = np.abs(E_std[:, 1])
    x_spread = np.abs(E_std[:, 4])
    x_switch = np.abs(E_std[:, 6])
    X = np.column_stack([np.ones(len(z)), x_std, x_spread, x_switch])
    beta, *_ = np.linalg.lstsq(X, z, rcond=None)
    return beta


def predict_sigma(beta: np.ndarray, E_std: np.ndarray) -> np.ndarray:
    x_std = np.abs(E_std[:, 1])
    x_spread = np.abs(E_std[:, 4])
    x_switch = np.abs(E_std[:, 6])
    X = np.column_stack([np.ones(E_std.shape[0]), x_std, x_spread, x_switch])
    log_var = X @ beta
    sigma = np.sqrt(np.exp(np.clip(log_var, -20.0, 8.0)))
    return np.clip(sigma, 1e-6, 10.0)


def hetero_loglike(y: np.ndarray, yhat: np.ndarray, sigma: np.ndarray) -> float:
    z = (y - yhat) / sigma
    return float(-0.5 * np.sum(z * z + 2.0 * np.log(sigma)))


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
        raise RuntimeError(f"Heteroscedastic OOS cohort too small: {len(rows)}")

    theta = np.array([r["theta_deg"] for r in rows], dtype=float)
    y = np.array([r["hxy"] for r in rows], dtype=float)
    hd = np.clip(helper.hellings_downs(theta), 1e-9, None)

    feature_names = ["f_mean", "f_std", "f_slope", "f_quad", "f_spread", "f_autoc1", "f_switch"]
    E_raw = np.column_stack([np.array([r[k] for r in rows], dtype=float) for k in feature_names])

    q_grid = np.linspace(max(0.9, q_center - 0.25), min(2.3, q_center + 0.25), 17)
    alpha_grid = np.array([0.08, 0.15, 0.30, 0.60, 1.20, 2.40], dtype=float)

    rng = np.random.default_rng(19701)
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

        # Baseline M2
        q2, beta2, sigma2 = m1817.best_q_m2(hd_tr, y_tr, q_grid=q_grid)
        X2_te = np.column_stack([hd_te ** q2, np.ones_like(hd_te)])
        pred2_te = X2_te @ beta2
        ll2_te = m1817.gaussian_loglike(y_te, pred2_te, sigma2)

        # Mean M2E model (as in QW-1817)
        qE, alphaE, betaE, sigmaE = m1817.best_q_alpha_m2e(hd_tr, y_tr, E_tr, q_grid=q_grid, alpha_grid=alpha_grid, rng=rng)
        XE_tr = np.column_stack([hd_tr ** qE, E_tr, np.ones_like(hd_tr)])
        XE_te = np.column_stack([hd_te ** qE, E_te, np.ones_like(hd_te)])
        predE_tr = XE_tr @ betaE
        predE_te = XE_te @ betaE

        # Homoscedastic LL (reference)
        llE_homo_te = m1817.gaussian_loglike(y_te, predE_te, sigmaE)

        # Heteroscedastic LL using variance model conditioned on embedding
        beta_var = fit_log_variance_model(y_tr - predE_tr, E_tr)
        sigma_tr_vec = predict_sigma(beta_var, E_tr)
        sigma_te_vec = predict_sigma(beta_var, E_te)

        # Calibrate overall scale on train to avoid systematic under/overconfidence
        resid_tr = y_tr - predE_tr
        scale = np.sqrt(np.mean((resid_tr / sigma_tr_vec) ** 2))
        scale = max(float(scale), 1e-6)
        sigma_te_cal = sigma_te_vec * scale

        llE_hetero_te = hetero_loglike(y_te, predE_te, sigma_te_cal)

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
                "ll_test_m2e_homo": float(llE_homo_te),
                "ll_test_m2e_hetero": float(llE_hetero_te),
                "delta_ll_hetero_vs_homo": float(llE_hetero_te - llE_homo_te),
                "delta_ll_hetero_vs_m2": float(llE_hetero_te - ll2_te),
                "delta_ll_homo_vs_m2": float(llE_homo_te - ll2_te),
            }
        )

    if len(rep_rows) < 10:
        raise RuntimeError("Too few heteroscedastic OOS replications.")

    arr_homo_vs_m2 = np.array([r["delta_ll_homo_vs_m2"] for r in rep_rows], dtype=float)
    arr_het_vs_m2 = np.array([r["delta_ll_hetero_vs_m2"] for r in rep_rows], dtype=float)
    arr_het_vs_homo = np.array([r["delta_ll_hetero_vs_homo"] for r in rep_rows], dtype=float)

    std_prev = float(d1817["summary"]["std_delta_ll_m2e_vs_m2"])
    std_homo_now = float(np.std(arr_homo_vs_m2))
    std_het_now = float(np.std(arr_het_vs_m2))

    summary = {
        "n_pairs_oos_cohort": len(rows),
        "mean_windows_per_pair": float(np.mean([r["n_windows"] for r in rows])),
        "n_features": int(E_raw.shape[1]),
        "n_replications": len(rep_rows),
        "train_fraction": 0.72,
        "mean_delta_ll_hetero_vs_homo": float(np.mean(arr_het_vs_homo)),
        "prob_hetero_better_than_homo": float(np.mean(arr_het_vs_homo > 0.0)),
        "mean_delta_ll_hetero_vs_m2": float(np.mean(arr_het_vs_m2)),
        "prob_hetero_better_than_m2": float(np.mean(arr_het_vs_m2 > 0.0)),
        "std_delta_ll_homo_vs_m2": std_homo_now,
        "std_delta_ll_hetero_vs_m2": std_het_now,
        "std_prev_qw1817": std_prev,
        "std_reduction_vs_qw1817": float(std_prev - std_het_now),
        "std_reduction_ratio_vs_qw1817": float((std_prev - std_het_now) / max(std_prev, 1e-12)),
        "std_reduction_hetero_vs_homo": float(std_homo_now - std_het_now),
    }

    pass_gain = summary["mean_delta_ll_hetero_vs_homo"] > 0.0 and summary["prob_hetero_better_than_homo"] >= 0.75
    pass_vs_m2 = summary["mean_delta_ll_hetero_vs_m2"] > 0.0 and summary["prob_hetero_better_than_m2"] >= 0.80
    pass_disp = summary["std_reduction_ratio_vs_qw1817"] >= 0.20

    if pass_gain and pass_vs_m2 and pass_disp:
        verdict = "HETEROSCEDASTIC_SEQUENCE_OOS_SUPPORTED"
    elif pass_gain and pass_vs_m2:
        verdict = "HETEROSCEDASTIC_SEQUENCE_OOS_PARTIAL"
    else:
        verdict = "HETEROSCEDASTIC_SEQUENCE_OOS_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "feature_names": feature_names,
        "summary": summary,
        "pass_flags": {
            "hetero_gain_vs_homo": bool(pass_gain),
            "hetero_gain_vs_m2": bool(pass_vs_m2),
            "dispersion_reduction_vs_qw1817": bool(pass_disp),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1820: HETEROSCEDASTIC SEQUENCE OOS",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- OOS cohort size: {len(rows)}",
        f"- Mean windows per pair: {summary['mean_windows_per_pair']:.2f}",
        f"- Features: {', '.join(feature_names)}",
        f"- Replications: {summary['n_replications']}",
        f"- Mean delta LL (hetero-homo): {summary['mean_delta_ll_hetero_vs_homo']:.4f}",
        f"- P(hetero>homo): {summary['prob_hetero_better_than_homo']:.3f}",
        f"- Mean delta LL (hetero-M2): {summary['mean_delta_ll_hetero_vs_m2']:.4f}",
        f"- P(hetero>M2): {summary['prob_hetero_better_than_m2']:.3f}",
        f"- Std delta LL vs M2 (QW-1817 -> hetero): {summary['std_prev_qw1817']:.3f} -> {summary['std_delta_ll_hetero_vs_m2']:.3f}",
        f"- Std reduction ratio: {summary['std_reduction_ratio_vs_qw1817']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- hetero_gain_vs_homo: {pass_gain}",
        f"- hetero_gain_vs_m2: {pass_vs_m2}",
        f"- dispersion_reduction_vs_qw1817: {pass_disp}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1820] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1820] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
