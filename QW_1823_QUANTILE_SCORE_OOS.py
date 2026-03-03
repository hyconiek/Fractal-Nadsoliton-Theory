#!/usr/bin/env python3
"""
QW-1823: Quantile-score OOS validation (distribution-free).

Implements robust predictive scoring requested after QW-1821/1822:
- no Gaussian/Student-t likelihood assumptions,
- uses empirical residual quantiles on train,
- evaluates pinball loss on test at taus = [0.1, 0.5, 0.9].

Lower score is better; we report gain = score(M2) - score(M2E).
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1823_quantile_score_oos.json"
OUT_MD = ROOT / "RAPORT_QW1823_QUANTILE_SCORE_OOS.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pinball_loss(y: np.ndarray, qhat: np.ndarray, tau: float) -> float:
    e = y - qhat
    return float(np.mean(np.where(e >= 0.0, tau * e, (tau - 1.0) * e)))


def quantile_score(y: np.ndarray, mu: np.ndarray, resid_train: np.ndarray, taus: np.ndarray) -> float:
    qs = []
    for tau in taus:
        q_res = float(np.quantile(resid_train, tau))
        qhat = mu + q_res
        qs.append(pinball_loss(y, qhat, float(tau)))
    return float(np.mean(qs))


def main() -> None:
    helper = load_module(ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py", "qw1787_helper")
    m1817 = load_module(ROOT / "QW_1817_SEQUENCE_OOS_VALIDATION.py", "qw1817_oos")

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
        raise RuntimeError(f"Quantile OOS cohort too small: {len(rows)}")

    theta = np.array([r["theta_deg"] for r in rows], dtype=float)
    y = np.array([r["hxy"] for r in rows], dtype=float)
    hd = np.clip(helper.hellings_downs(theta), 1e-9, None)

    feature_names = ["f_mean", "f_std", "f_slope", "f_quad", "f_spread", "f_autoc1", "f_switch"]
    E_raw = np.column_stack([np.array([r[k] for r in rows], dtype=float) for k in feature_names])

    q_grid = np.linspace(max(0.9, q_center - 0.25), min(2.3, q_center + 0.25), 17)
    alpha_grid = np.array([0.08, 0.15, 0.30, 0.60, 1.20, 2.40], dtype=float)
    taus = np.array([0.10, 0.50, 0.90], dtype=float)

    rng = np.random.default_rng(19901)
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

        # M2 model
        q2, beta2, _ = m1817.best_q_m2(hd_tr, y_tr, q_grid=q_grid)
        X2_tr = np.column_stack([hd_tr ** q2, np.ones_like(hd_tr)])
        X2_te = np.column_stack([hd_te ** q2, np.ones_like(hd_te)])
        mu2_tr = X2_tr @ beta2
        mu2_te = X2_te @ beta2
        r2_tr = y_tr - mu2_tr

        # M2E model
        qE, alphaE, betaE, _ = m1817.best_q_alpha_m2e(hd_tr, y_tr, E_tr, q_grid=q_grid, alpha_grid=alpha_grid, rng=rng)
        XE_tr = np.column_stack([hd_tr ** qE, E_tr, np.ones_like(hd_tr)])
        XE_te = np.column_stack([hd_te ** qE, E_te, np.ones_like(hd_te)])
        muE_tr = XE_tr @ betaE
        muE_te = XE_te @ betaE
        rE_tr = y_tr - muE_tr

        qs2 = quantile_score(y_te, mu2_te, r2_tr, taus=taus)
        qsE = quantile_score(y_te, muE_te, rE_tr, taus=taus)

        mae2 = float(np.mean(np.abs(y_te - mu2_te)))
        maeE = float(np.mean(np.abs(y_te - muE_te)))

        rep_rows.append(
            {
                "rep": i,
                "n_train": int(len(tr_idx)),
                "n_test": int(len(te_idx)),
                "q_m2": float(q2),
                "q_m2e": float(qE),
                "alpha_m2e": float(alphaE),
                "quantile_score_m2": float(qs2),
                "quantile_score_m2e": float(qsE),
                "quantile_gain_m2_minus_m2e": float(qs2 - qsE),
                "mae_m2": mae2,
                "mae_m2e": maeE,
                "mae_gain_m2_minus_m2e": float(mae2 - maeE),
            }
        )

    if len(rep_rows) < 10:
        raise RuntimeError("Too few quantile-score OOS replications.")

    arr_q_gain = np.array([r["quantile_gain_m2_minus_m2e"] for r in rep_rows], dtype=float)
    arr_mae_gain = np.array([r["mae_gain_m2_minus_m2e"] for r in rep_rows], dtype=float)

    summary = {
        "n_pairs_oos_cohort": len(rows),
        "mean_windows_per_pair": float(np.mean([r["n_windows"] for r in rows])),
        "n_features": int(E_raw.shape[1]),
        "n_replications": len(rep_rows),
        "train_fraction": 0.72,
        "taus": taus.tolist(),
        "mean_quantile_gain_m2_minus_m2e": float(np.mean(arr_q_gain)),
        "median_quantile_gain_m2_minus_m2e": float(np.median(arr_q_gain)),
        "std_quantile_gain_m2_minus_m2e": float(np.std(arr_q_gain)),
        "prob_quantile_gain_positive": float(np.mean(arr_q_gain > 0.0)),
        "mean_mae_gain_m2_minus_m2e": float(np.mean(arr_mae_gain)),
        "prob_mae_gain_positive": float(np.mean(arr_mae_gain > 0.0)),
    }

    pass_q = summary["mean_quantile_gain_m2_minus_m2e"] > 0.0 and summary["prob_quantile_gain_positive"] >= 0.80
    pass_mae = summary["mean_mae_gain_m2_minus_m2e"] > 0.0 and summary["prob_mae_gain_positive"] >= 0.90
    pass_disp = summary["std_quantile_gain_m2_minus_m2e"] <= 0.040

    if pass_q and pass_mae and pass_disp:
        verdict = "QUANTILE_SCORE_OOS_SUPPORTED"
    elif pass_q and pass_mae:
        verdict = "QUANTILE_SCORE_OOS_PARTIAL"
    else:
        verdict = "QUANTILE_SCORE_OOS_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "feature_names": feature_names,
        "summary": summary,
        "pass_flags": {
            "quantile_gain": bool(pass_q),
            "mae_gain": bool(pass_mae),
            "dispersion_control": bool(pass_disp),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1823: QUANTILE SCORE OOS",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- OOS cohort size: {len(rows)}",
        f"- Mean windows per pair: {summary['mean_windows_per_pair']:.2f}",
        f"- Features: {', '.join(feature_names)}",
        f"- Taus: {summary['taus']}",
        f"- Replications: {summary['n_replications']}",
        f"- Mean quantile gain (M2-M2E): {summary['mean_quantile_gain_m2_minus_m2e']:.6f}",
        f"- P(quantile gain>0): {summary['prob_quantile_gain_positive']:.3f}",
        f"- Mean MAE gain (M2-M2E): {summary['mean_mae_gain_m2_minus_m2e']:.6f}",
        f"- P(MAE gain>0): {summary['prob_mae_gain_positive']:.3f}",
        f"- Std quantile gain: {summary['std_quantile_gain_m2_minus_m2e']:.6f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- quantile_gain: {pass_q}",
        f"- mae_gain: {pass_mae}",
        f"- dispersion_control: {pass_disp}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1823] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1823] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
