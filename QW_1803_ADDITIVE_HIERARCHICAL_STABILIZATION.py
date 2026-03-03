#!/usr/bin/env python3
"""
QW-1803: Hierarchical stabilization for additive decoherence model.

Follow-up to QW-1802:
- additive form gives strong gain and better residual flattening,
- but replication dispersion is too high.
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1803_additive_hierarchical_stabilization.json"
OUT_MD = ROOT / "RAPORT_QW1803_ADDITIVE_HIERARCHICAL_STABILIZATION.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def damping_quadratic(theta_deg: np.ndarray, lam: float) -> np.ndarray:
    x = (1.0 - np.cos(np.radians(theta_deg))) / 2.0
    return np.exp(-lam * (x * x))


def additive_feature(theta_deg: np.ndarray, lam: float) -> np.ndarray:
    d = damping_quadratic(theta_deg, lam)
    g = 1.0 - d
    return g - float(np.mean(g))


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


def evidence_m6_add(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int, *,
                    b_mean: float | None = None, b_std: float | None = None,
                    lam_mean: float | None = None, lam_std: float | None = None) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    if b_mean is None or b_std is None:
        B = rng.uniform(-1.2, 1.2, n_mc)
    else:
        B = np.clip(rng.normal(loc=b_mean, scale=max(b_std, 1e-3), size=n_mc), -1.2, 1.2)
    if lam_mean is None or lam_std is None:
        lam = rng.uniform(0.0, 6.0, n_mc)
    else:
        lam = np.clip(rng.normal(loc=lam_mean, scale=max(lam_std, 1e-3), size=n_mc), 0.0, 6.0)

    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    cache: Dict[float, np.ndarray] = {}
    ll_rows = []
    for a, b, c, qq, ll in zip(A, B, C, q, lam):
        key = float(np.round(ll, 3))
        if key not in cache:
            cache[key] = additive_feature(theta, key)
        g = cache[key]
        ll_rows.append(loglike(H, sigma, a * (hd0 ** qq) + b * g + c))
    ll = np.array(ll_rows, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def posterior_hyper_additive(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.2, 1.2, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    lam = rng.uniform(0.0, 6.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    cache: Dict[float, np.ndarray] = {}
    ll_rows = []
    for a, b, c, qq, ll in zip(A, B, C, q, lam):
        key = float(np.round(ll, 3))
        if key not in cache:
            cache[key] = additive_feature(theta, key)
        g = cache[key]
        ll_rows.append(loglike(H, sigma, a * (hd0 ** qq) + b * g + c))
    ll = np.array(ll_rows, dtype=float)
    lw = ll - logsumexp(ll)
    w = np.exp(lw)

    b_mean = float(np.sum(w * B))
    l_mean = float(np.sum(w * lam))
    b_var = float(np.sum(w * (B - b_mean) ** 2))
    l_var = float(np.sum(w * (lam - l_mean) ** 2))
    return {
        "b_mean": b_mean,
        "lam_mean": l_mean,
        "b_std_raw": float(np.sqrt(max(b_var, 0.0))),
        "lam_std_raw": float(np.sqrt(max(l_var, 0.0))),
    }


def main() -> None:
    helper = load_module(ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py", "qw1787_helper")

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
        if len(x) >= n_match_min and stab <= stability_max:
            rows.append({"theta_deg": float(sep), "hxy": float(hxy)})

    if len(rows) < 85:
        raise RuntimeError("Operational cohort too small.")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)

    # Baseline full.
    z0 = evidence_flat(H_all, n_mc=8500, seed=18301)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=18302)
    full_m2 = float(z2 - z0)

    raw = posterior_hyper_additive(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=65000, seed=18303)
    b_std_raw = raw["b_std_raw"]
    lam_std_raw = raw["lam_std_raw"]

    shrink_grid = [0.20, 0.30, 0.40, 0.55, 0.75]
    sweep = []
    for j, sf in enumerate(shrink_grid):
        b_std = max(0.04, sf * b_std_raw)
        lam_std = max(0.08, sf * lam_std_raw)
        z6h = evidence_m6_add(
            helper,
            theta_all,
            H_all,
            q_center=q_center,
            q_width=q_width,
            n_mc=17000,
            seed=18400 + 100 * j + 1,
            b_mean=raw["b_mean"],
            b_std=b_std,
            lam_mean=raw["lam_mean"],
            lam_std=lam_std,
        )
        full_m6h = float(z6h - z0)
        full_delta = float(full_m6h - full_m2)

        rng = np.random.default_rng(18400 + 100 * j + 2)
        rep_rows = []
        for i in range(12):
            idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
            if len(idx) < 80:
                continue
            th = theta_all[idx]
            hh = H_all[idx]
            b0 = evidence_flat(hh, n_mc=3900, seed=18400 + 100 * j + 20 * i + 3)
            b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=18400 + 100 * j + 20 * i + 4)
            b6 = evidence_m6_add(
                helper,
                th,
                hh,
                q_center=q_center,
                q_width=q_width,
                n_mc=8200,
                seed=18400 + 100 * j + 20 * i + 5,
                b_mean=raw["b_mean"],
                b_std=b_std,
                lam_mean=raw["lam_mean"],
                lam_std=lam_std,
            )
            l2 = float(b2 - b0)
            l6 = float(b6 - b0)
            rep_rows.append({"delta_hier_vs_m2": float(l6 - l2), "logB_hier_vs_flat": l6})

        arr_d = np.array([r["delta_hier_vs_m2"] for r in rep_rows], dtype=float)
        arr_h = np.array([r["logB_hier_vs_flat"] for r in rep_rows], dtype=float)
        p_gain = float(np.mean(arr_d > 0.0))
        p_flat = float(np.mean(arr_h > 0.0))
        std_d = float(np.std(arr_d))
        med_d = float(np.median(arr_d))

        pass_basic = full_delta > 0.0 and p_gain >= 0.90 and p_flat >= 0.95 and std_d <= 0.35
        score = (
            0.35 * min(1.0, full_delta / 1.5)
            + 0.25 * p_gain
            + 0.20 * p_flat
            + 0.20 * max(0.0, 1.0 - std_d / 0.50)
        )
        sweep.append(
            {
                "shrink_factor": sf,
                "b_std": b_std,
                "lam_std": lam_std,
                "full_logB_m2_vs_flat": full_m2,
                "full_logB_hier_vs_flat": full_m6h,
                "full_delta_hier_vs_m2": full_delta,
                "replications": len(rep_rows),
                "prob_hier_gt_m2": p_gain,
                "prob_hier_gt_flat": p_flat,
                "median_delta_hier_vs_m2": med_d,
                "std_delta_hier_vs_m2": std_d,
                "pass_basic": bool(pass_basic),
                "selection_score": float(score),
            }
        )

    valid = [r for r in sweep if r["pass_basic"]]
    if len(valid) > 0:
        best = max(valid, key=lambda r: r["selection_score"])
        recommendation_strength = "STRONG"
    else:
        best = max(sweep, key=lambda r: r["selection_score"])
        recommendation_strength = "CONDITIONAL"

    if recommendation_strength == "STRONG":
        verdict = "ADDITIVE_HIERARCHICAL_STABILIZATION_SUPPORTED"
    else:
        verdict = "ADDITIVE_HIERARCHICAL_STABILIZATION_PARTIAL"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_pairs": len(rows),
        "fraction": frac,
        "q_width": q_width,
        "raw_posterior": raw,
        "sweep": sweep,
        "best": best,
        "recommendation_strength": recommendation_strength,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1803: ADDITIVE HIERARCHICAL STABILIZATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Best shrink factor: {best['shrink_factor']:.2f} ({recommendation_strength})",
        f"- Full delta hier-M2 (best): {best['full_delta_hier_vs_m2']:.4f}",
        f"- P(hier>M2) (best): {best['prob_hier_gt_m2']:.3f}",
        f"- P(hier>flat) (best): {best['prob_hier_gt_flat']:.3f}",
        f"- Std delta hier-M2 (best): {best['std_delta_hier_vs_m2']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Sweep",
    ]
    for r in sweep:
        lines.append(
            f"- sf={r['shrink_factor']:.2f}: full_delta={r['full_delta_hier_vs_m2']:.4f}, P(hier>M2)={r['prob_hier_gt_m2']:.3f}, "
            f"P(hier>flat)={r['prob_hier_gt_flat']:.3f}, std_delta={r['std_delta_hier_vs_m2']:.3f}, "
            f"score={r['selection_score']:.3f}, pass_basic={r['pass_basic']}"
        )
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1803] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1803] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
