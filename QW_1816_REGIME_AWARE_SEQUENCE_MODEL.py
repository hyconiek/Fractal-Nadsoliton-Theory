#!/usr/bin/env python3
"""
QW-1816: Regime-aware sequence model on top of multiscale embedding.

Goal:
- keep the strong gain from QW-1815,
- explicitly model latent heterogeneity across pairs via regime dummies,
- test whether replication dispersion can be reduced.

Model:
  M2ER(theta, E, R) = A * HD(theta)^q + <B, E> + g1*R1 + g2*R2 + C
where R1,R2 are regime indicators derived from PC1 tertiles of embedding E.
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


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1816_regime_aware_sequence_model.json"
OUT_MD = ROOT / "RAPORT_QW1816_REGIME_AWARE_SEQUENCE_MODEL.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def zscore(v: np.ndarray) -> np.ndarray:
    s = float(np.std(v))
    if s <= 1e-12:
        return np.zeros_like(v)
    return (v - float(np.mean(v))) / s


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


def build_regime_dummies(E: np.ndarray) -> Tuple[np.ndarray, Dict[str, object]]:
    """
    Build 3 balanced latent regimes from PC1(E), then return two dummies:
    R1 = regime 1, R2 = regime 2 (regime 0 is baseline).
    """
    Ec = E - np.mean(E, axis=0, keepdims=True)
    _, _, vt = np.linalg.svd(Ec, full_matrices=False)
    pc1 = Ec @ vt[0]

    q1, q2 = np.quantile(pc1, [1.0 / 3.0, 2.0 / 3.0])
    reg = np.zeros(len(pc1), dtype=int)
    reg[pc1 > q1] = 1
    reg[pc1 > q2] = 2

    r1 = (reg == 1).astype(float)
    r2 = (reg == 2).astype(float)
    counts = [int(np.sum(reg == k)) for k in (0, 1, 2)]

    meta = {
        "pc1_tertiles": [float(q1), float(q2)],
        "regime_counts": counts,
        "regime_fractions": [float(c / len(reg)) for c in counts],
    }
    return np.column_stack([r1, r2]), meta


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


def evidence_m2e(helper, theta: np.ndarray, H: np.ndarray, E: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    k = E.shape[1]
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    tau = rng.uniform(0.05, 0.35, n_mc)
    B = np.array([rng.normal(0.0, t, size=k) for t in tau], dtype=float)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)

    ll_rows = []
    for i in range(n_mc):
        base = A[i] * (hd0 ** q[i]) + C[i]
        add = E @ B[i]
        ll_rows.append(loglike(H, sigma, base + add))
    ll = np.array(ll_rows, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_m2er(
    helper,
    theta: np.ndarray,
    H: np.ndarray,
    E: np.ndarray,
    R: np.ndarray,
    q_center: float,
    q_width: float,
    n_mc: int,
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    k = E.shape[1]
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)

    tau_b = rng.uniform(0.05, 0.30, n_mc)
    B = np.array([rng.normal(0.0, t, size=k) for t in tau_b], dtype=float)

    tau_g = rng.uniform(0.02, 0.40, n_mc)
    G = np.array([rng.normal(0.0, t, size=2) for t in tau_g], dtype=float)

    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)

    ll_rows = []
    for i in range(n_mc):
        base = A[i] * (hd0 ** q[i]) + C[i]
        add = E @ B[i]
        reg = R @ G[i]
        ll_rows.append(loglike(H, sigma, base + add + reg))
    ll = np.array(ll_rows, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


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
        raise RuntimeError(f"Regime-aware cohort too small: {len(rows)}")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    feature_names = ["f_mean", "f_std", "f_slope", "f_quad", "f_spread", "f_autoc1", "f_switch"]
    E_raw = np.column_stack([np.array([r[k] for r in rows], dtype=float) for k in feature_names])
    E = np.column_stack([zscore(E_raw[:, i]) for i in range(E_raw.shape[1])])

    R, regime_meta = build_regime_dummies(E)

    z0 = evidence_flat(H_all, n_mc=8500, seed=19401)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=19402)
    ze = evidence_m2e(helper, theta_all, H_all, E, q_center=q_center, q_width=q_width, n_mc=15000, seed=19403)
    zer = evidence_m2er(helper, theta_all, H_all, E, R, q_center=q_center, q_width=q_width, n_mc=18000, seed=19404)

    full_m2 = float(z2 - z0)
    full_m2e = float(ze - z0)
    full_m2er = float(zer - z0)

    full_delta_er_vs_m2 = float(full_m2er - full_m2)
    full_delta_er_vs_m2e = float(full_m2er - full_m2e)

    rng = np.random.default_rng(19405)
    rep_rows = []
    for i in range(12):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 75:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        ee = E[idx]
        rr = R[idx]

        b0 = evidence_flat(hh, n_mc=3900, seed=19410 + 20 * i + 1)
        b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=19410 + 20 * i + 2)
        be = evidence_m2e(helper, th, hh, ee, q_center=q_center, q_width=q_width, n_mc=7600, seed=19410 + 20 * i + 3)
        ber = evidence_m2er(helper, th, hh, ee, rr, q_center=q_center, q_width=q_width, n_mc=9800, seed=19410 + 20 * i + 4)

        l2 = float(b2 - b0)
        le = float(be - b0)
        ler = float(ber - b0)
        rep_rows.append(
            {
                "rep": i,
                "logB_m2_vs_flat": l2,
                "logB_m2e_vs_flat": le,
                "logB_m2er_vs_flat": ler,
                "delta_m2e_vs_m2": float(le - l2),
                "delta_m2er_vs_m2": float(ler - l2),
                "delta_m2er_vs_m2e": float(ler - le),
            }
        )

    arr_d_e = np.array([r["delta_m2e_vs_m2"] for r in rep_rows], dtype=float)
    arr_d_er = np.array([r["delta_m2er_vs_m2"] for r in rep_rows], dtype=float)
    arr_er_vs_e = np.array([r["delta_m2er_vs_m2e"] for r in rep_rows], dtype=float)
    arr_er = np.array([r["logB_m2er_vs_flat"] for r in rep_rows], dtype=float)

    std_e = float(np.std(arr_d_e))
    std_er = float(np.std(arr_d_er))
    std_reduction = float(std_e - std_er)

    summary = {
        "n_pairs_regime_cohort": len(rows),
        "mean_windows_per_pair": float(np.mean([r["n_windows"] for r in rows])),
        "n_features": E.shape[1],
        "regime_meta": regime_meta,
        "fraction": frac,
        "q_width": q_width,
        "full_logB_m2_vs_flat": full_m2,
        "full_logB_m2e_vs_flat": full_m2e,
        "full_logB_m2er_vs_flat": full_m2er,
        "full_delta_m2er_vs_m2": full_delta_er_vs_m2,
        "full_delta_m2er_vs_m2e": full_delta_er_vs_m2e,
        "replications": len(rep_rows),
        "prob_m2er_gt_m2": float(np.mean(arr_d_er > 0.0)),
        "prob_m2er_gt_m2e": float(np.mean(arr_er_vs_e > 0.0)),
        "prob_m2er_gt_flat": float(np.mean(arr_er > 0.0)),
        "median_delta_m2er_vs_m2": float(np.median(arr_d_er)),
        "std_delta_m2e_vs_m2": std_e,
        "std_delta_m2er_vs_m2": std_er,
        "std_reduction_vs_m2e": std_reduction,
    }

    pass_full = summary["full_delta_m2er_vs_m2e"] > 0.0
    pass_rep = summary["prob_m2er_gt_m2e"] >= 0.80 and summary["prob_m2er_gt_flat"] >= 0.95
    pass_disp = summary["std_delta_m2er_vs_m2"] <= 0.30
    pass_stab = summary["std_reduction_vs_m2e"] >= 0.10

    if pass_full and pass_rep and pass_disp and pass_stab:
        verdict = "REGIME_AWARE_SEQUENCE_MODEL_SUPPORTED"
    elif pass_full and pass_rep and pass_stab:
        verdict = "REGIME_AWARE_SEQUENCE_MODEL_PARTIAL"
    else:
        verdict = "REGIME_AWARE_SEQUENCE_MODEL_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "feature_names": feature_names,
        "summary": summary,
        "pass_flags": {
            "full_gain_vs_m2e": bool(pass_full),
            "replication_gain_vs_m2e": bool(pass_rep),
            "dispersion_control": bool(pass_disp),
            "stability_improvement_vs_m2e": bool(pass_stab),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1816: REGIME-AWARE SEQUENCE MODEL",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Regime-aware cohort size: {len(rows)}",
        f"- Mean windows per pair: {summary['mean_windows_per_pair']:.2f}",
        f"- Features: {', '.join(feature_names)}",
        f"- Regime counts: {regime_meta['regime_counts']}",
        f"- Full logB M2/M2E/M2ER: {full_m2:.4f} / {full_m2e:.4f} / {full_m2er:.4f}",
        f"- Full delta M2ER-M2E: {full_delta_er_vs_m2e:.4f}",
        f"- P(M2ER>M2E): {summary['prob_m2er_gt_m2e']:.3f}",
        f"- P(M2ER>flat): {summary['prob_m2er_gt_flat']:.3f}",
        f"- Std delta (M2E-M2): {std_e:.3f}",
        f"- Std delta (M2ER-M2): {std_er:.3f}",
        f"- Std reduction vs M2E: {std_reduction:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_gain_vs_m2e: {pass_full}",
        f"- replication_gain_vs_m2e: {pass_rep}",
        f"- dispersion_control: {pass_disp}",
        f"- stability_improvement_vs_m2e: {pass_stab}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1816] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1816] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
