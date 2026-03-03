#!/usr/bin/env python3
"""
QW-1848: PTA unit-of-analysis audit for quantile-gain probability criterion.

Compares:
- split-level positivity (used in QW-1823 summary),
- pair-level positivity under repeated OOS cross-splits.

This detects potential pseudo-replication / averaging compression effects.
"""

from __future__ import annotations

import importlib.util
import json
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import beta, binom


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1848_pta_unit_of_analysis_audit.json"
OUT_MD = ROOT / "RAPORT_QW1848_PTA_UNIT_OF_ANALYSIS_AUDIT.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pinball_single(y: float, qhat: float, tau: float) -> float:
    e = y - qhat
    return float(tau * e if e >= 0.0 else (tau - 1.0) * e)


def quantile_score_single(y: float, mu: float, resid_train: np.ndarray, taus: np.ndarray) -> float:
    vals = []
    for tau in taus:
        q_res = float(np.quantile(resid_train, tau))
        vals.append(pinball_single(y, mu + q_res, float(tau)))
    return float(np.mean(vals))


def binom_tail(n: int, k: int, p0: float) -> float:
    return float(binom.sf(k - 1, n, p0))


def one_sided_lower_bound(k: int, n: int, alpha: float = 0.05) -> float:
    if k <= 0:
        return 0.0
    # Exact Clopper-Pearson lower one-sided bound
    return float(beta.ppf(alpha, k, n - k + 1))


def build_rows(helper, m1817, n_match_min: int, stability_max: float) -> List[Dict[str, float]]:
    residuals = helper.load_residuals(
        ROOT / "nano15/residuals/NANOGrav15yr_PulsarTiming_v2.1.0/residuals",
        max_psr=34,
    )
    positions = helper.load_positions(ROOT / "nano15/parfiles")

    rows: List[Dict[str, float]] = []
    for p1, p2 in combinations(list(residuals.keys()), 2):
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
            }
        )
    return rows


def main() -> None:
    helper = load_module(ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py", "qw1787_helper")
    m1817 = load_module(ROOT / "QW_1817_SEQUENCE_OOS_VALIDATION.py", "qw1817_oos")

    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))

    q_center = float(d1773["projection"]["p"])
    cohort = d1793["operational_protocol"]["cohort"]
    n_match_min = int(cohort["n_match_min"])
    stability_max = float(cohort["stability_max"])

    rows = build_rows(helper, m1817, n_match_min=n_match_min, stability_max=stability_max)
    if len(rows) < 80:
        raise RuntimeError(f"PTA cohort too small for audit: {len(rows)}")

    theta = np.array([r["theta_deg"] for r in rows], dtype=float)
    y = np.array([r["hxy"] for r in rows], dtype=float)
    hd = np.clip(helper.hellings_downs(theta), 1e-9, None)

    feature_names = ["f_mean", "f_std", "f_slope", "f_quad", "f_spread", "f_autoc1", "f_switch"]
    E_raw = np.column_stack([np.array([r[k] for r in rows], dtype=float) for k in feature_names])

    q_grid = np.linspace(max(0.9, q_center - 0.25), min(2.3, q_center + 0.25), 17)
    alpha_grid = np.array([0.08, 0.15, 0.30, 0.60, 1.20, 2.40], dtype=float)
    taus = np.array([0.10, 0.50, 0.90], dtype=float)

    rng = np.random.default_rng(21930)
    n_rep = 80

    pair_gains: List[List[float]] = [[] for _ in range(len(rows))]
    rep_rows: List[Dict[str, float]] = []

    for rep in range(n_rep):
        tr_idx, te_idx = m1817.stratified_split_indices(theta, train_frac=0.72, rng=rng, n_bins=8)
        if len(tr_idx) < 55 or len(te_idx) < 20:
            continue

        hd_tr, hd_te = hd[tr_idx], hd[te_idx]
        y_tr, y_te = y[tr_idx], y[te_idx]
        E_tr_raw, E_te_raw = E_raw[tr_idx], E_raw[te_idx]
        E_tr, E_te = m1817.standardize_train_test(E_tr_raw, E_te_raw)

        q2, beta2, _ = m1817.best_q_m2(hd_tr, y_tr, q_grid=q_grid)
        X2_tr = np.column_stack([hd_tr ** q2, np.ones_like(hd_tr)])
        X2_te = np.column_stack([hd_te ** q2, np.ones_like(hd_te)])
        mu2_tr = X2_tr @ beta2
        mu2_te = X2_te @ beta2
        r2_tr = y_tr - mu2_tr

        qE, alphaE, betaE, _ = m1817.best_q_alpha_m2e(
            hd_tr,
            y_tr,
            E_tr,
            q_grid=q_grid,
            alpha_grid=alpha_grid,
            rng=np.random.default_rng(50000 + rep),
        )
        XE_tr = np.column_stack([hd_tr ** qE, E_tr, np.ones_like(hd_tr)])
        XE_te = np.column_stack([hd_te ** qE, E_te, np.ones_like(hd_te)])
        muE_tr = XE_tr @ betaE
        muE_te = XE_te @ betaE
        rE_tr = y_tr - muE_tr

        gains = []
        for j, idx in enumerate(te_idx):
            qs2 = quantile_score_single(float(y_te[j]), float(mu2_te[j]), r2_tr, taus=taus)
            qsE = quantile_score_single(float(y_te[j]), float(muE_te[j]), rE_tr, taus=taus)
            g = float(qs2 - qsE)
            gains.append(g)
            pair_gains[int(idx)].append(g)

        gains_arr = np.array(gains, dtype=float)
        rep_rows.append(
            {
                "rep": int(rep),
                "n_test": int(len(te_idx)),
                "q_m2": float(q2),
                "q_m2e": float(qE),
                "alpha_m2e": float(alphaE),
                "mean_pair_gain": float(np.mean(gains_arr)),
                "prob_pair_gain_positive": float(np.mean(gains_arr > 0.0)),
            }
        )

    if len(rep_rows) < 30:
        raise RuntimeError("Too few audit replications computed.")

    pair_obs = np.array([len(v) for v in pair_gains if len(v) > 0], dtype=int)
    pair_mean_gain = np.array([float(np.mean(v)) for v in pair_gains if len(v) > 0], dtype=float)
    all_unit_gain = np.concatenate([np.array(v, dtype=float) for v in pair_gains if len(v) > 0])

    rep_mean_gain = np.array([r["mean_pair_gain"] for r in rep_rows], dtype=float)
    rep_prob_pos = np.array([r["prob_pair_gain_positive"] for r in rep_rows], dtype=float)

    k_pair_mean_pos = int(np.sum(pair_mean_gain > 0.0))
    n_pair = int(len(pair_mean_gain))

    pval_pair_vs_0p9 = binom_tail(n_pair, k_pair_mean_pos, 0.90)
    pval_pair_vs_0p8 = binom_tail(n_pair, k_pair_mean_pos, 0.80)
    lower95_pair = one_sided_lower_bound(k_pair_mean_pos, n_pair, alpha=0.05)

    split_level_prob = float(np.mean(rep_mean_gain > 0.0))
    pair_level_prob = float(k_pair_mean_pos / n_pair)

    summary = {
        "n_pairs_cohort": int(len(rows)),
        "n_pairs_observed": n_pair,
        "n_replications": int(len(rep_rows)),
        "pair_obs_per_index": {
            "min": int(np.min(pair_obs)),
            "median": float(np.median(pair_obs)),
            "max": int(np.max(pair_obs)),
        },
        "split_level": {
            "prob_rep_mean_gain_positive": split_level_prob,
            "mean_rep_mean_gain": float(np.mean(rep_mean_gain)),
            "std_rep_mean_gain": float(np.std(rep_mean_gain)),
            "mean_rep_prob_pair_positive": float(np.mean(rep_prob_pos)),
        },
        "pair_level": {
            "prob_pair_mean_gain_positive": pair_level_prob,
            "mean_pair_mean_gain": float(np.mean(pair_mean_gain)),
            "std_pair_mean_gain": float(np.std(pair_mean_gain)),
            "prob_all_unit_gain_positive": float(np.mean(all_unit_gain > 0.0)),
            "mean_all_unit_gain": float(np.mean(all_unit_gain)),
            "std_all_unit_gain": float(np.std(all_unit_gain)),
        },
        "inference": {
            "k_pair_mean_gain_positive": k_pair_mean_pos,
            "n_pairs_tested": n_pair,
            "lower95_one_sided_for_prob_pair_mean_positive": lower95_pair,
            "pvalue_h0_prob_le_0p90": pval_pair_vs_0p9,
            "pvalue_h0_prob_le_0p80": pval_pair_vs_0p8,
        },
        "compression_gap": float(split_level_prob - pair_level_prob),
    }

    pass_split = split_level_prob >= 0.90
    pass_pair_0p9 = lower95_pair >= 0.90 and pval_pair_vs_0p9 <= 0.05
    pass_pair_0p8 = lower95_pair >= 0.80 and pval_pair_vs_0p8 <= 0.05

    if pass_split and (not pass_pair_0p9):
        verdict = "PTA_UNIT_MISMATCH_REQUIRES_CRITERION_REDESIGN"
    elif pass_pair_0p9:
        verdict = "PTA_UNIT_CONSISTENT_STRONG"
    elif pass_pair_0p8:
        verdict = "PTA_UNIT_CONSISTENT_MODERATE"
    else:
        verdict = "PTA_UNIT_CONSISTENT_WEAK"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "pass_flags": {
            "split_level_prob_ge_0p9": bool(pass_split),
            "pair_level_prob_ge_0p9_inferential": bool(pass_pair_0p9),
            "pair_level_prob_ge_0p8_inferential": bool(pass_pair_0p8),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    s = summary
    lines = [
        "# RAPORT QW-1848: PTA UNIT OF ANALYSIS AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Cohort pairs / observed pairs: {s['n_pairs_cohort']} / {s['n_pairs_observed']}",
        f"- Replications: {s['n_replications']}",
        (
            "- Pair obs per index (min/med/max): "
            f"{s['pair_obs_per_index']['min']} / {s['pair_obs_per_index']['median']:.1f} / {s['pair_obs_per_index']['max']}"
        ),
        "",
        "## Split vs Pair",
        f"- Split-level P(rep_mean_gain>0): {s['split_level']['prob_rep_mean_gain_positive']:.3f}",
        f"- Pair-level P(pair_mean_gain>0): {s['pair_level']['prob_pair_mean_gain_positive']:.3f}",
        f"- Compression gap: {s['compression_gap']:.3f}",
        "",
        "## Pair-Level Inference",
        (
            f"- k/n positives: {s['inference']['k_pair_mean_gain_positive']}/"
            f"{s['inference']['n_pairs_tested']}"
        ),
        (
            "- one-sided lower95 for prob: "
            f"{s['inference']['lower95_one_sided_for_prob_pair_mean_positive']:.3f}"
        ),
        f"- p-value H0: prob<=0.90 -> {s['inference']['pvalue_h0_prob_le_0p90']:.6f}",
        f"- p-value H0: prob<=0.80 -> {s['inference']['pvalue_h0_prob_le_0p80']:.6f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- split_level_prob_ge_0p9: {pass_split}",
        f"- pair_level_prob_ge_0p9_inferential: {pass_pair_0p9}",
        f"- pair_level_prob_ge_0p8_inferential: {pass_pair_0p8}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1848] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1848] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
