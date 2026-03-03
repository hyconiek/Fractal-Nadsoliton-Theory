#!/usr/bin/env python3
"""
QW-1797: Hierarchical multimode shrinkage test.

Purpose:
- follow up QW-1796 where multimode extension had local gain but high variance,
- use shared (hierarchical) priors for extension parameters (B, alpha),
- test if replication stability improves without adding new model complexity.
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
from numpy.polynomial import legendre as npleg
from scipy.special import logsumexp
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1797_hierarchical_multimode_shrinkage.json"
OUT_MD = ROOT / "RAPORT_QW1797_HIERARCHICAL_MULTIMODE_SHRINKAGE.md"


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
    weights = np.array([(ell + 1.0) ** (-alpha) for ell in ell_list], dtype=float)
    weights /= np.sum(np.abs(weights))
    T = np.zeros_like(theta_deg, dtype=float)
    for w, ell in zip(weights, ell_list):
        T += w * legendre(theta_deg, ell)
    s = float(np.std(T))
    if s > 1e-12:
        T /= s
    return T


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


def evidence_m5(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, ell_list: List[int], n_mc: int, seed: int, *,
                b_prior_mean: float | None = None, b_prior_std: float | None = None,
                alpha_prior_mean: float | None = None, alpha_prior_std: float | None = None) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)

    if b_prior_mean is None or b_prior_std is None:
        B = rng.uniform(-1.0, 1.0, n_mc)
    else:
        B = np.clip(rng.normal(loc=b_prior_mean, scale=max(b_prior_std, 1e-3), size=n_mc), -1.0, 1.0)

    if alpha_prior_mean is None or alpha_prior_std is None:
        alpha = rng.uniform(0.5, 3.0, n_mc)
    else:
        alpha = np.clip(rng.normal(loc=alpha_prior_mean, scale=max(alpha_prior_std, 1e-3), size=n_mc), 0.5, 3.0)

    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)

    tmpl_cache: Dict[float, np.ndarray] = {}
    ll_rows = []
    for a, b, c, qq, aa in zip(A, B, C, q, alpha):
        key = float(np.round(aa, 3))
        if key not in tmpl_cache:
            tmpl_cache[key] = template_family(theta, ell_list, key)
        T = tmpl_cache[key]
        ll_rows.append(loglike(H, sigma, a * (hd0 ** qq) + b * T + c))

    ll = np.array(ll_rows, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def posterior_hyper_from_full(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, ell_list: List[int], n_mc: int, seed: int) -> Dict[str, float]:
    """
    Approximate posterior moments for B and alpha from full cohort under broad priors.
    """
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
        ll_rows.append(loglike(H, sigma, a * (hd0 ** qq) + b * T + c))
    ll = np.array(ll_rows, dtype=float)

    # normalize weights in log-space
    lw = ll - logsumexp(ll)
    w = np.exp(lw)

    b_mean = float(np.sum(w * B))
    a_mean = float(np.sum(w * alpha))
    b_var = float(np.sum(w * (B - b_mean) ** 2))
    a_var = float(np.sum(w * (alpha - a_mean) ** 2))

    # shrink priors: half of posterior std, with floors
    b_std = max(0.05, 0.5 * np.sqrt(max(b_var, 1e-8)))
    a_std = max(0.08, 0.5 * np.sqrt(max(a_var, 1e-8)))
    return {
        "b_mean": b_mean,
        "b_std": b_std,
        "alpha_mean": a_mean,
        "alpha_std": a_std,
        "posterior_b_std_raw": float(np.sqrt(max(b_var, 0.0))),
        "posterior_alpha_std_raw": float(np.sqrt(max(a_var, 0.0))),
    }


def fit_best_m2(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    best_ll = -np.inf
    best = None
    for a, c, qq in zip(A, C, q):
        m = a * (hd0 ** qq) + c
        ll = loglike(H, sigma, m)
        if ll > best_ll:
            best_ll = ll
            best = {"A": float(a), "C": float(c), "q": float(qq), "sigma": sigma, "loglike": best_ll}
    return best


def fit_best_m5_hier(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, ell_list: List[int], prior: Dict[str, float], n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    B = np.clip(rng.normal(loc=prior["b_mean"], scale=prior["b_std"], size=n_mc), -1.0, 1.0)
    alpha = np.clip(rng.normal(loc=prior["alpha_mean"], scale=prior["alpha_std"], size=n_mc), 0.5, 3.0)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)

    tmpl_cache: Dict[float, np.ndarray] = {}
    best_ll = -np.inf
    best = None
    for a, b, c, qq, aa in zip(A, B, C, q, alpha):
        key = float(np.round(aa, 3))
        if key not in tmpl_cache:
            tmpl_cache[key] = template_family(theta, ell_list, key)
        T = tmpl_cache[key]
        m = a * (hd0 ** qq) + b * T + c
        ll = loglike(H, sigma, m)
        if ll > best_ll:
            best_ll = ll
            best = {"A": float(a), "B": float(b), "C": float(c), "q": float(qq), "alpha": float(aa), "sigma": sigma, "loglike": best_ll}
    return best


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


def main() -> None:
    helper = load_helper()
    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))
    d1796 = json.loads((ROOT / "report_qw1796_physical_multimode_extension.json").read_text(encoding="utf-8"))

    q_center = float(d1773["projection"]["p"])
    q_width = float(d1793["operational_protocol"]["q_width"])
    frac = float(d1793["operational_protocol"]["fraction"])
    cohort = d1793["operational_protocol"]["cohort"]
    n_match_min = int(cohort["n_match_min"])
    stability_max = float(cohort["stability_max"])

    # Continue with best family from QW-1796.
    best_family = d1796["best_family"]
    family_name = str(best_family["family"])
    ell_list = list(best_family["ell_list"])

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

    # Learn hierarchical priors for extension parameters from full cohort.
    prior = posterior_hyper_from_full(
        helper, theta_all, H_all, q_center=q_center, q_width=q_width, ell_list=ell_list, n_mc=70000, seed=17971
    )

    # Full-cohort evidence
    z0 = evidence_flat(H_all, n_mc=8500, seed=17972)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=17973)
    z5_local = evidence_m5(
        helper, theta_all, H_all, q_center=q_center, q_width=q_width, ell_list=ell_list, n_mc=17000, seed=17974
    )
    z5_hier = evidence_m5(
        helper,
        theta_all,
        H_all,
        q_center=q_center,
        q_width=q_width,
        ell_list=ell_list,
        n_mc=17000,
        seed=17975,
        b_prior_mean=prior["b_mean"],
        b_prior_std=prior["b_std"],
        alpha_prior_mean=prior["alpha_mean"],
        alpha_prior_std=prior["alpha_std"],
    )

    full_m2 = float(z2 - z0)
    full_m5_local = float(z5_local - z0)
    full_m5_hier = float(z5_hier - z0)
    full_delta_hier_vs_m2 = float(full_m5_hier - full_m2)
    full_delta_hier_vs_local = float(full_m5_hier - full_m5_local)

    # Residual flattening
    fit2 = fit_best_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=80000, seed=17976)
    fit5h = fit_best_m5_hier(
        helper, theta_all, H_all, q_center=q_center, q_width=q_width, ell_list=ell_list, prior=prior, n_mc=100000, seed=17977
    )
    hd0 = np.clip(helper.hellings_downs(theta_all), 1e-9, None)
    pred2 = fit2["A"] * (hd0 ** fit2["q"]) + fit2["C"]
    T5h = template_family(theta_all, ell_list, fit5h["alpha"])
    pred5h = fit5h["A"] * (hd0 ** fit5h["q"]) + fit5h["B"] * T5h + fit5h["C"]
    rs2 = residual_stats(H_all - pred2, theta_all)
    rs5h = residual_stats(H_all - pred5h, theta_all)

    # Replications
    rng = np.random.default_rng(17978)
    rep_rows = []
    for i in range(14):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 80:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        b0 = evidence_flat(hh, n_mc=4200, seed=17980 + 30 * i + 1)
        b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6800, seed=17980 + 30 * i + 2)
        b5h = evidence_m5(
            helper,
            th,
            hh,
            q_center=q_center,
            q_width=q_width,
            ell_list=ell_list,
            n_mc=9000,
            seed=17980 + 30 * i + 3,
            b_prior_mean=prior["b_mean"],
            b_prior_std=prior["b_std"],
            alpha_prior_mean=prior["alpha_mean"],
            alpha_prior_std=prior["alpha_std"],
        )
        l2 = float(b2 - b0)
        l5h = float(b5h - b0)
        rep_rows.append({"rep": i, "logB_m2_vs_flat": l2, "logB_m5_hier_vs_flat": l5h, "delta_hier_vs_m2": float(l5h - l2)})

    arr_h = np.array([r["logB_m5_hier_vs_flat"] for r in rep_rows], dtype=float)
    arr_d = np.array([r["delta_hier_vs_m2"] for r in rep_rows], dtype=float)

    summary = {
        "family": family_name,
        "ell_list": ell_list,
        "n_pairs": len(rows),
        "fraction": frac,
        "q_width": q_width,
        "full_logB_m2_vs_flat": full_m2,
        "full_logB_m5_local_vs_flat": full_m5_local,
        "full_logB_m5_hier_vs_flat": full_m5_hier,
        "full_delta_hier_vs_m2": full_delta_hier_vs_m2,
        "full_delta_hier_vs_local": full_delta_hier_vs_local,
        "replications": len(rep_rows),
        "prob_hier_gt_flat": float(np.mean(arr_h > 0.0)),
        "prob_hier_gt_m2": float(np.mean(arr_d > 0.0)),
        "median_delta_hier_vs_m2": float(np.median(arr_d)),
        "std_delta_hier_vs_m2": float(np.std(arr_d)),
        "resid_m2_rho_theta": rs2["rho_theta_resid"],
        "resid_hier_rho_theta": rs5h["rho_theta_resid"],
        "resid_m2_max_abs_bin_mean": rs2["max_abs_bin_mean"],
        "resid_hier_max_abs_bin_mean": rs5h["max_abs_bin_mean"],
        "resid_rho_improvement_abs": float(abs(rs2["rho_theta_resid"]) - abs(rs5h["rho_theta_resid"])),
        "resid_bin_improvement": float(rs2["max_abs_bin_mean"] - rs5h["max_abs_bin_mean"]),
    }

    pass_full = summary["full_delta_hier_vs_m2"] > 0.0
    pass_rep = summary["prob_hier_gt_m2"] >= 0.70 and summary["prob_hier_gt_flat"] >= 0.90
    pass_disp = summary["std_delta_hier_vs_m2"] <= 0.22
    pass_resid = summary["resid_rho_improvement_abs"] > 0.02 or summary["resid_bin_improvement"] > 0.03

    if pass_full and pass_rep and pass_disp and pass_resid:
        verdict = "HIERARCHICAL_MULTIMODE_SUPPORTED"
    elif pass_full and pass_rep:
        verdict = "HIERARCHICAL_MULTIMODE_PARTIAL"
    else:
        verdict = "HIERARCHICAL_MULTIMODE_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "hierarchical_prior": prior,
        "pass_flags": {
            "full_gain_over_m2": bool(pass_full),
            "replication_gain": bool(pass_rep),
            "dispersion_control": bool(pass_disp),
            "residual_improvement": bool(pass_resid),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1797: HIERARCHICAL MULTIMODE SHRINKAGE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Family: {family_name} {ell_list}",
        f"- Full logB M2/local/hier: {full_m2:.4f} / {full_m5_local:.4f} / {full_m5_hier:.4f}",
        f"- Full delta hier-M2: {full_delta_hier_vs_m2:.4f}",
        f"- P(hier>M2): {summary['prob_hier_gt_m2']:.3f}",
        f"- P(hier>flat): {summary['prob_hier_gt_flat']:.3f}",
        f"- Std delta hier-M2: {summary['std_delta_hier_vs_m2']:.3f}",
        f"- Residual improvements: d|rho|={summary['resid_rho_improvement_abs']:.4f}, d(bin)={summary['resid_bin_improvement']:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Hierarchical Prior",
        f"- b_mean={prior['b_mean']:.4f}, b_std={prior['b_std']:.4f} (raw_post_std={prior['posterior_b_std_raw']:.4f})",
        f"- alpha_mean={prior['alpha_mean']:.4f}, alpha_std={prior['alpha_std']:.4f} (raw_post_std={prior['posterior_alpha_std_raw']:.4f})",
        "",
        "## Pass Flags",
        f"- full_gain_over_m2: {pass_full}",
        f"- replication_gain: {pass_rep}",
        f"- dispersion_control: {pass_disp}",
        f"- residual_improvement: {pass_resid}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1797] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1797] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
