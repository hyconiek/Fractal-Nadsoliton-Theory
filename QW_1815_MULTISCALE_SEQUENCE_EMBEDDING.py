#!/usr/bin/env python3
"""
QW-1815: Multiscale sequence embedding model.

Redesign after QW-1814:
- move from scalar dynamic proxies to richer multiscale sequence descriptors.

Model:
  M2E(theta, E) = A * HD(theta)^q + <B, E> + C
where E is a standardized embedding extracted from windowed cross-DFA trajectories.
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
OUT_JSON = ROOT / "report_qw1815_multiscale_sequence_embedding.json"
OUT_MD = ROOT / "RAPORT_QW1815_MULTISCALE_SEQUENCE_EMBEDDING.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def center(v: np.ndarray) -> np.ndarray:
    return v - float(np.mean(v))


def zscore(v: np.ndarray) -> np.ndarray:
    s = float(np.std(v))
    if s <= 1e-12:
        return np.zeros_like(v)
    return (v - float(np.mean(v))) / s


def windowed_multiscale_traj(helper, x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns:
      centers: normalized window centers in [0,1]
      vals:    cross-DFA values for windows
    """
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

    # linear + quadratic trend in window-center time
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
    # shrinked priors on embedding coefficients for stability
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
        raise RuntimeError(f"Embedding cohort too small: {len(rows)}")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    feature_names = ["f_mean", "f_std", "f_slope", "f_quad", "f_spread", "f_autoc1", "f_switch"]
    E_raw = np.column_stack([np.array([r[k] for r in rows], dtype=float) for k in feature_names])
    E = np.column_stack([zscore(E_raw[:, i]) for i in range(E_raw.shape[1])])

    z0 = evidence_flat(H_all, n_mc=8500, seed=19301)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=19302)
    ze = evidence_m2e(helper, theta_all, H_all, E, q_center=q_center, q_width=q_width, n_mc=17000, seed=19303)
    full_m2 = float(z2 - z0)
    full_m2e = float(ze - z0)
    full_delta = float(full_m2e - full_m2)

    rng = np.random.default_rng(19304)
    rep_rows = []
    for i in range(12):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 75:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        ee = E[idx]
        b0 = evidence_flat(hh, n_mc=3900, seed=19310 + 20 * i + 1)
        b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=19310 + 20 * i + 2)
        be = evidence_m2e(helper, th, hh, ee, q_center=q_center, q_width=q_width, n_mc=8600, seed=19310 + 20 * i + 3)
        l2 = float(b2 - b0)
        le = float(be - b0)
        rep_rows.append({"rep": i, "logB_m2_vs_flat": l2, "logB_m2e_vs_flat": le, "delta_m2e_vs_m2": float(le - l2)})

    arr_d = np.array([r["delta_m2e_vs_m2"] for r in rep_rows], dtype=float)
    arr_e = np.array([r["logB_m2e_vs_flat"] for r in rep_rows], dtype=float)

    summary = {
        "n_pairs_embedding_cohort": len(rows),
        "mean_windows_per_pair": float(np.mean([r["n_windows"] for r in rows])),
        "n_features": E.shape[1],
        "fraction": frac,
        "q_width": q_width,
        "full_logB_m2_vs_flat": full_m2,
        "full_logB_m2e_vs_flat": full_m2e,
        "full_delta_m2e_vs_m2": full_delta,
        "replications": len(rep_rows),
        "prob_m2e_gt_m2": float(np.mean(arr_d > 0.0)),
        "prob_m2e_gt_flat": float(np.mean(arr_e > 0.0)),
        "median_delta_m2e_vs_m2": float(np.median(arr_d)),
        "std_delta_m2e_vs_m2": float(np.std(arr_d)),
    }

    pass_full = summary["full_delta_m2e_vs_m2"] > 0.0
    pass_rep = summary["prob_m2e_gt_m2"] >= 0.80 and summary["prob_m2e_gt_flat"] >= 0.95
    pass_disp = summary["std_delta_m2e_vs_m2"] <= 0.30

    if pass_full and pass_rep and pass_disp:
        verdict = "MULTISCALE_SEQUENCE_EMBEDDING_SUPPORTED"
    elif pass_full and pass_rep:
        verdict = "MULTISCALE_SEQUENCE_EMBEDDING_PARTIAL"
    else:
        verdict = "MULTISCALE_SEQUENCE_EMBEDDING_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "feature_names": feature_names,
        "summary": summary,
        "pass_flags": {
            "full_gain": bool(pass_full),
            "replication_gain": bool(pass_rep),
            "dispersion_control": bool(pass_disp),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1815: MULTISCALE SEQUENCE EMBEDDING",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Embedding cohort size: {len(rows)}",
        f"- Mean windows per pair: {summary['mean_windows_per_pair']:.2f}",
        f"- Features: {', '.join(feature_names)}",
        f"- Full logB M2/M2E: {full_m2:.4f} / {full_m2e:.4f}",
        f"- Full delta M2E-M2: {full_delta:.4f}",
        f"- P(M2E>M2): {summary['prob_m2e_gt_m2']:.3f}",
        f"- P(M2E>flat): {summary['prob_m2e_gt_flat']:.3f}",
        f"- Std delta M2E-M2: {summary['std_delta_m2e_vs_m2']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_gain: {pass_full}",
        f"- replication_gain: {pass_rep}",
        f"- dispersion_control: {pass_disp}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1815] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1815] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
