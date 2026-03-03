#!/usr/bin/env python3
"""
QW-1796: Physically constrained multimode angular extension.

Baseline (locked by QW-1793):
    M2(theta) = A * HD(theta)^q + C

Extension families (2 new dof only: B, alpha):
    M5(theta) = A * HD(theta)^q + B * T_family(theta; alpha) + C

where
    T_family(theta; alpha) = sum_l w_l(alpha) P_l(cos theta),
    w_l(alpha) ~ (l+1)^(-alpha), alpha in [0.5, 3.0].

Goal:
- check if a physically constrained multimode component can improve evidence
  and reduce residual angular structure without overfitting.
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from numpy.polynomial import legendre as npleg
from scipy.stats import spearmanr
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1796_physical_multimode_extension.json"
OUT_MD = ROOT / "RAPORT_QW1796_PHYSICAL_MULTIMODE_EXTENSION.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def legendre(theta_deg: np.ndarray, ell: int) -> np.ndarray:
    x = np.cos(np.radians(theta_deg))
    coeff = np.zeros(ell + 1, dtype=float)
    coeff[ell] = 1.0
    return npleg.legval(x, coeff)


def template_family(theta_deg: np.ndarray, ell_list: List[int], alpha: float) -> np.ndarray:
    terms = []
    weights = []
    for ell in ell_list:
        terms.append(legendre(theta_deg, ell))
        weights.append((ell + 1.0) ** (-alpha))
    W = np.array(weights, dtype=float)
    W /= np.sum(np.abs(W))
    T = np.zeros_like(theta_deg, dtype=float)
    for w, t in zip(W, terms):
        T += w * t
    # normalize to unit std for stable B interpretation
    s = float(np.std(T))
    if s > 1e-12:
        T /= s
    return T


def loglike_gauss(H: np.ndarray, sigma: float, model: np.ndarray) -> float:
    z = (H - model) / sigma
    return float(-0.5 * np.sum(z * z))


def fit_best_m2(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)

    best = None
    best_ll = -np.inf
    for a, c, qq in zip(A, C, q):
        m = a * (hd0 ** qq) + c
        ll = loglike_gauss(H, sigma, m)
        if ll > best_ll:
            best_ll = ll
            best = {"A": float(a), "C": float(c), "q": float(qq), "sigma": sigma, "loglike": best_ll}
    return best


def fit_best_m5(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, ell_list: List[int], n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.0, 1.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    alpha = rng.uniform(0.5, 3.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)

    # cache templates by rounded alpha for speed
    tmpl_cache: Dict[float, np.ndarray] = {}
    best = None
    best_ll = -np.inf
    for a, b, c, qq, aa in zip(A, B, C, q, alpha):
        key = float(np.round(aa, 3))
        if key not in tmpl_cache:
            tmpl_cache[key] = template_family(theta, ell_list, key)
        T = tmpl_cache[key]
        m = a * (hd0 ** qq) + b * T + c
        ll = loglike_gauss(H, sigma, m)
        if ll > best_ll:
            best_ll = ll
            best = {
                "A": float(a),
                "B": float(b),
                "C": float(c),
                "q": float(qq),
                "alpha": float(aa),
                "sigma": sigma,
                "loglike": best_ll,
            }
    return best


def evidence_m2(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array([loglike_gauss(H, sigma, a * (hd0 ** qq) + c) for a, c, qq in zip(A, C, q)], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_m5(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, ell_list: List[int], n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.0, 1.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    alpha = rng.uniform(0.5, 3.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)

    tmpl_cache: Dict[float, np.ndarray] = {}
    ll_rows = []
    for a, b, c, qq, aa in zip(A, B, C, q, alpha):
        key = float(np.round(aa, 3))
        if key not in tmpl_cache:
            tmpl_cache[key] = template_family(theta, ell_list, key)
        T = tmpl_cache[key]
        ll_rows.append(loglike_gauss(H, sigma, a * (hd0 ** qq) + b * T + c))
    ll = np.array(ll_rows, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_flat(H: np.ndarray, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 2.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    ll = np.array([loglike_gauss(H, sigma, np.full_like(H, c)) for c in C], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def residual_stats(resid: np.ndarray, theta: np.ndarray) -> Dict[str, float]:
    rho, _ = spearmanr(theta, resid)
    rho = float(rho)
    bins = np.linspace(0.0, 180.0, 7)
    means = []
    for i in range(len(bins) - 1):
        m = (theta >= bins[i]) & (theta < bins[i + 1] if i < len(bins) - 2 else theta <= bins[i + 1])
        if np.sum(m) == 0:
            continue
        means.append(float(np.mean(resid[m])))
    max_bin = float(max(abs(v) for v in means)) if means else 0.0
    return {"rho_theta_resid": rho, "max_abs_bin_mean": max_bin}


def select_best_family(results: List[Dict[str, object]]) -> Dict[str, object]:
    # prefer families with positive full gain and stronger replication support
    return max(
        results,
        key=lambda r: (
            float(r["full_delta_m5_minus_m2"]),
            float(r["prob_m5_gt_m2"]),
            -float(r["std_delta_m5_minus_m2"]),
        ),
    )


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
        if len(x) >= n_match_min and stab <= stability_max:
            rows.append({"theta_deg": float(sep), "hxy": float(hxy)})

    if len(rows) < 85:
        raise RuntimeError("Operational cohort too small.")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)

    families = [
        {"name": "odd_multiscale", "ell_list": [1, 3, 5]},
        {"name": "even_multiscale", "ell_list": [2, 4, 6]},
        {"name": "mixed_low", "ell_list": [1, 2, 3]},
    ]

    family_results: List[Dict[str, object]] = []
    for j, fam in enumerate(families):
        name = fam["name"]
        ell_list = fam["ell_list"]

        z0 = evidence_flat(H_all, n_mc=8000, seed=19600 + 100 * j + 1)
        z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=12000, seed=19600 + 100 * j + 2)
        z5 = evidence_m5(helper, theta_all, H_all, q_center=q_center, q_width=q_width, ell_list=ell_list, n_mc=18000, seed=19600 + 100 * j + 3)

        full_m2 = float(z2 - z0)
        full_m5 = float(z5 - z0)
        full_delta = float(full_m5 - full_m2)

        # Residual flattening (full data, best fits).
        fit2 = fit_best_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=90000, seed=19600 + 100 * j + 4)
        fit5 = fit_best_m5(helper, theta_all, H_all, q_center=q_center, q_width=q_width, ell_list=ell_list, n_mc=110000, seed=19600 + 100 * j + 5)

        hd0 = np.clip(helper.hellings_downs(theta_all), 1e-9, None)
        pred2 = fit2["A"] * (hd0 ** fit2["q"]) + fit2["C"]
        T5 = template_family(theta_all, ell_list, fit5["alpha"])
        pred5 = fit5["A"] * (hd0 ** fit5["q"]) + fit5["B"] * T5 + fit5["C"]
        rs2 = residual_stats(H_all - pred2, theta_all)
        rs5 = residual_stats(H_all - pred5, theta_all)

        # Replications
        rng = np.random.default_rng(17960 + j)
        rep_rows = []
        for i in range(12):
            idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
            if len(idx) < 80:
                continue
            th = theta_all[idx]
            hh = H_all[idx]
            b0 = evidence_flat(hh, n_mc=4500, seed=19700 + 100 * j + 20 * i + 1)
            b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=7000, seed=19700 + 100 * j + 20 * i + 2)
            b5 = evidence_m5(helper, th, hh, q_center=q_center, q_width=q_width, ell_list=ell_list, n_mc=9500, seed=19700 + 100 * j + 20 * i + 3)
            l2 = float(b2 - b0)
            l5 = float(b5 - b0)
            rep_rows.append({"rep": i, "logB_m2_vs_flat": l2, "logB_m5_vs_flat": l5, "delta_m5_minus_m2": float(l5 - l2)})

        arr_m5 = np.array([r["logB_m5_vs_flat"] for r in rep_rows], dtype=float)
        arr_d = np.array([r["delta_m5_minus_m2"] for r in rep_rows], dtype=float)

        family_results.append(
            {
                "family": name,
                "ell_list": ell_list,
                "full_logB_m2_vs_flat": full_m2,
                "full_logB_m5_vs_flat": full_m5,
                "full_delta_m5_minus_m2": full_delta,
                "replications": len(rep_rows),
                "prob_m5_gt_flat": float(np.mean(arr_m5 > 0.0)),
                "prob_m5_gt_m2": float(np.mean(arr_d > 0.0)),
                "median_delta_m5_minus_m2": float(np.median(arr_d)),
                "std_delta_m5_minus_m2": float(np.std(arr_d)),
                "resid_m2_rho_theta": rs2["rho_theta_resid"],
                "resid_m5_rho_theta": rs5["rho_theta_resid"],
                "resid_m2_max_abs_bin_mean": rs2["max_abs_bin_mean"],
                "resid_m5_max_abs_bin_mean": rs5["max_abs_bin_mean"],
                "resid_rho_improvement_abs": float(abs(rs2["rho_theta_resid"]) - abs(rs5["rho_theta_resid"])),
                "resid_bin_improvement": float(rs2["max_abs_bin_mean"] - rs5["max_abs_bin_mean"]),
            }
        )

    best = select_best_family(family_results)

    pass_full = float(best["full_delta_m5_minus_m2"]) > 0.0
    pass_rep = float(best["prob_m5_gt_m2"]) >= 0.75 and float(best["prob_m5_gt_flat"]) >= 0.95
    pass_disp = float(best["std_delta_m5_minus_m2"]) <= 0.18
    pass_resid = float(best["resid_rho_improvement_abs"]) > 0.02 and float(best["resid_bin_improvement"]) > 0.03

    if pass_full and pass_rep and pass_disp and pass_resid:
        verdict = "PHYSICAL_MULTIMODE_EXTENSION_SUPPORTED"
    elif pass_full and pass_rep:
        verdict = "PHYSICAL_MULTIMODE_EXTENSION_PARTIAL"
    else:
        verdict = "PHYSICAL_MULTIMODE_EXTENSION_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "locked_protocol_used": {
            "fraction": frac,
            "q_width": q_width,
            "cohort_name": cohort.get("name"),
            "n_pairs": len(rows),
        },
        "family_results": family_results,
        "best_family": best,
        "pass_flags": {
            "full_gain_over_m2": bool(pass_full),
            "replication_gain": bool(pass_rep),
            "dispersion_control": bool(pass_disp),
            "residual_flattening": bool(pass_resid),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1796: PHYSICAL MULTIMODE EXTENSION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Locked protocol: frac={frac}, q_width={q_width}, cohort={cohort.get('name')}, n={len(rows)}",
        f"- Best family: {best['family']} (ells={best['ell_list']})",
        f"- Full delta (M5-M2): {best['full_delta_m5_minus_m2']:.4f}",
        f"- P(M5>M2): {best['prob_m5_gt_m2']:.3f}",
        f"- P(M5>flat): {best['prob_m5_gt_flat']:.3f}",
        f"- Residual improvements: d|rho|={best['resid_rho_improvement_abs']:.4f}, d(bin)={best['resid_bin_improvement']:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_gain_over_m2: {pass_full}",
        f"- replication_gain: {pass_rep}",
        f"- dispersion_control: {pass_disp}",
        f"- residual_flattening: {pass_resid}",
        "",
        "## Family Table",
    ]
    for r in family_results:
        lines.append(
            f"- {r['family']} {r['ell_list']}: full_delta={r['full_delta_m5_minus_m2']:.4f}, "
            f"P(M5>M2)={r['prob_m5_gt_m2']:.3f}, P(M5>flat)={r['prob_m5_gt_flat']:.3f}, "
            f"std_delta={r['std_delta_m5_minus_m2']:.3f}, d|rho|={r['resid_rho_improvement_abs']:.4f}, "
            f"d(bin)={r['resid_bin_improvement']:.4f}"
        )
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1796] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1796] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
