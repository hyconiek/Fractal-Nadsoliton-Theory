#!/usr/bin/env python3
"""
QW-1809: Dynamic feature scan (single-feature latent models).

Scans three dynamic proxies:
  F1: drift      = h2 - h1
  F2: volatility = |h2 - h1|
  F3: log_ratio  = log((|h2|+eps)/(|h1|+eps))

Model family:
  M2F(theta, F) = A * HD(theta)^q + B * F + C
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
OUT_JSON = ROOT / "report_qw1809_dynamic_feature_scan.json"
OUT_MD = ROOT / "RAPORT_QW1809_DYNAMIC_FEATURE_SCAN.md"


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


def evidence_m2f(helper, theta: np.ndarray, H: np.ndarray, F: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.5, 1.5, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array([loglike(H, sigma, a * (hd0 ** qq) + b * F + c) for a, b, c, qq in zip(A, B, C, q)], dtype=float)
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
        h1, h2 = split_half_hxy(helper, x, y)
        if not np.isfinite(h1) or not np.isfinite(h2):
            continue
        drift = float(h2 - h1)
        stability = float(abs(drift))
        if len(x) >= n_match_min and stability <= stability_max:
            eps = 1e-6
            rows.append(
                {
                    "theta_deg": float(sep),
                    "hxy": float(hxy),
                    "drift": drift,
                    "volatility": float(abs(drift)),
                    "log_ratio": float(np.log((abs(h2) + eps) / (abs(h1) + eps))),
                }
            )

    if len(rows) < 85:
        raise RuntimeError("Operational cohort too small.")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    feat_map = {
        "drift": center(np.array([r["drift"] for r in rows], dtype=float)),
        "volatility": center(np.array([r["volatility"] for r in rows], dtype=float)),
        "log_ratio": center(np.array([r["log_ratio"] for r in rows], dtype=float)),
    }

    z0 = evidence_flat(H_all, n_mc=8500, seed=18901)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=18902)
    full_m2 = float(z2 - z0)

    results = []
    for j, (fname, F_all) in enumerate(feat_map.items()):
        zf = evidence_m2f(helper, theta_all, H_all, F_all, q_center=q_center, q_width=q_width, n_mc=17000, seed=18910 + 100 * j + 1)
        full_m2f = float(zf - z0)
        full_delta = float(full_m2f - full_m2)

        rng = np.random.default_rng(18910 + 100 * j + 2)
        rep_rows = []
        for i in range(12):
            idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
            if len(idx) < 80:
                continue
            th = theta_all[idx]
            hh = H_all[idx]
            ff = F_all[idx]
            b0 = evidence_flat(hh, n_mc=3900, seed=18910 + 100 * j + 20 * i + 3)
            b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=18910 + 100 * j + 20 * i + 4)
            bf = evidence_m2f(helper, th, hh, ff, q_center=q_center, q_width=q_width, n_mc=8200, seed=18910 + 100 * j + 20 * i + 5)
            l2 = float(b2 - b0)
            lf = float(bf - b0)
            rep_rows.append({"delta_m2f_vs_m2": float(lf - l2), "logB_m2f_vs_flat": lf})

        arr_d = np.array([r["delta_m2f_vs_m2"] for r in rep_rows], dtype=float)
        arr_f = np.array([r["logB_m2f_vs_flat"] for r in rep_rows], dtype=float)
        results.append(
            {
                "feature": fname,
                "full_logB_m2_vs_flat": full_m2,
                "full_logB_m2f_vs_flat": full_m2f,
                "full_delta_m2f_vs_m2": full_delta,
                "replications": len(rep_rows),
                "prob_m2f_gt_m2": float(np.mean(arr_d > 0.0)),
                "prob_m2f_gt_flat": float(np.mean(arr_f > 0.0)),
                "median_delta_m2f_vs_m2": float(np.median(arr_d)),
                "std_delta_m2f_vs_m2": float(np.std(arr_d)),
            }
        )

    best = max(results, key=lambda r: (float(r["full_delta_m2f_vs_m2"]), float(r["prob_m2f_gt_m2"]), -float(r["std_delta_m2f_vs_m2"])))
    pass_gain = best["full_delta_m2f_vs_m2"] > 0.0 and best["prob_m2f_gt_m2"] >= 0.80 and best["prob_m2f_gt_flat"] >= 0.95
    pass_disp = best["std_delta_m2f_vs_m2"] <= 0.30

    if pass_gain and pass_disp:
        verdict = "DYNAMIC_FEATURE_SCAN_SUPPORTED"
    elif pass_gain:
        verdict = "DYNAMIC_FEATURE_SCAN_PARTIAL"
    else:
        verdict = "DYNAMIC_FEATURE_SCAN_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_pairs": len(rows),
        "fraction": frac,
        "q_width": q_width,
        "results": results,
        "best_feature": best,
        "pass_flags": {"gain": bool(pass_gain), "dispersion_control": bool(pass_disp)},
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1809: DYNAMIC FEATURE SCAN",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Best feature: {best['feature']}",
        f"- Full delta best: {best['full_delta_m2f_vs_m2']:.4f}",
        f"- P(best>M2): {best['prob_m2f_gt_m2']:.3f}",
        f"- P(best>flat): {best['prob_m2f_gt_flat']:.3f}",
        f"- Std delta best: {best['std_delta_m2f_vs_m2']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Results",
    ]
    for r in results:
        lines.append(
            f"- {r['feature']}: full_delta={r['full_delta_m2f_vs_m2']:.4f}, P(>M2)={r['prob_m2f_gt_m2']:.3f}, "
            f"P(>flat)={r['prob_m2f_gt_flat']:.3f}, std={r['std_delta_m2f_vs_m2']:.3f}"
        )
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1809] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1809] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
