#!/usr/bin/env python3
"""
QW-1801: Physical decoherence extension for PTA angular model.

Motivation:
- hierarchical angular-mode branch (QW-1796..1800) was not robust,
- residual structure remained angle-linked (QW-1794),
- test a different physical mechanism: propagation decoherence damping.

Baseline:
    M2(theta) = A * HD(theta)^q + C

Extension family:
    M6(theta) = A * HD(theta)^q * D(theta; lambda) + C

with one additional parameter lambda >= 0 and three constrained damping laws:
    linear:    D = exp(-lambda * x)
    quadratic: D = exp(-lambda * x^2)
    rational:  D = 1 / (1 + lambda * x)
where x = (1 - cos(theta))/2.
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
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1801_physical_decoherence_extension.json"
OUT_MD = ROOT / "RAPORT_QW1801_PHYSICAL_DECOHERENCE_EXTENSION.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_m2_fit_module():
    path = ROOT / "QW_1797_HIERARCHICAL_MULTIMODE_SHRINKAGE.py"
    spec = importlib.util.spec_from_file_location("qw1797_core", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def damping(theta_deg: np.ndarray, lam: float, kind: str) -> np.ndarray:
    x = (1.0 - np.cos(np.radians(theta_deg))) / 2.0
    if kind == "linear":
        return np.exp(-lam * x)
    if kind == "quadratic":
        return np.exp(-lam * (x * x))
    if kind == "rational":
        return 1.0 / (1.0 + lam * x)
    raise ValueError(f"Unknown damping kind: {kind}")


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


def evidence_m6(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, kind: str, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    lam = rng.uniform(0.0, 6.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)

    # Cache damping by rounded lambda to reduce recomputation.
    damp_cache: Dict[float, np.ndarray] = {}
    ll_rows = []
    for a, c, qq, ll in zip(A, C, q, lam):
        key = float(np.round(ll, 3))
        if key not in damp_cache:
            damp_cache[key] = damping(theta, key, kind)
        d = damp_cache[key]
        model = a * (hd0 ** qq) * d + c
        ll_rows.append(loglike(H, sigma, model))
    ll_arr = np.array(ll_rows, dtype=float)
    return float(logsumexp(ll_arr) - np.log(len(ll_arr)))


def fit_best_m6(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, kind: str, n_mc: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    lam = rng.uniform(0.0, 6.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)

    damp_cache: Dict[float, np.ndarray] = {}
    best_ll = -np.inf
    best = None
    for a, c, qq, ll in zip(A, C, q, lam):
        key = float(np.round(ll, 3))
        if key not in damp_cache:
            damp_cache[key] = damping(theta, key, kind)
        d = damp_cache[key]
        model = a * (hd0 ** qq) * d + c
        lk = loglike(H, sigma, model)
        if lk > best_ll:
            best_ll = lk
            best = {"A": float(a), "C": float(c), "q": float(qq), "lambda": float(ll), "sigma": sigma, "loglike": best_ll}
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
    m2mod = load_m2_fit_module()

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
    hd0 = np.clip(helper.hellings_downs(theta_all), 1e-9, None)

    # Baseline M2 on full cohort.
    z0 = evidence_flat(H_all, n_mc=8500, seed=18010)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=18011)
    full_m2 = float(z2 - z0)
    fit2 = m2mod.fit_best_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=70000, seed=18012)
    pred2 = fit2["A"] * (hd0 ** fit2["q"]) + fit2["C"]
    rs2 = residual_stats(H_all - pred2, theta_all)

    kinds = ["linear", "quadratic", "rational"]
    kind_results = []
    for j, kind in enumerate(kinds):
        z6 = evidence_m6(helper, theta_all, H_all, q_center=q_center, q_width=q_width, kind=kind, n_mc=18000, seed=18100 + 100 * j + 1)
        full_m6 = float(z6 - z0)
        full_delta = float(full_m6 - full_m2)

        fit6 = fit_best_m6(helper, theta_all, H_all, q_center=q_center, q_width=q_width, kind=kind, n_mc=95000, seed=18100 + 100 * j + 2)
        d6 = damping(theta_all, fit6["lambda"], kind)
        pred6 = fit6["A"] * (hd0 ** fit6["q"]) * d6 + fit6["C"]
        rs6 = residual_stats(H_all - pred6, theta_all)

        rng = np.random.default_rng(18100 + 100 * j + 3)
        rep_rows = []
        for i in range(12):
            idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
            if len(idx) < 80:
                continue
            th = theta_all[idx]
            hh = H_all[idx]
            b0 = evidence_flat(hh, n_mc=4200, seed=18200 + 100 * j + 20 * i + 1)
            b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6500, seed=18200 + 100 * j + 20 * i + 2)
            b6 = evidence_m6(helper, th, hh, q_center=q_center, q_width=q_width, kind=kind, n_mc=8500, seed=18200 + 100 * j + 20 * i + 3)
            l2 = float(b2 - b0)
            l6 = float(b6 - b0)
            rep_rows.append({"rep": i, "logB_m2_vs_flat": l2, "logB_m6_vs_flat": l6, "delta_m6_minus_m2": float(l6 - l2)})

        arr6 = np.array([r["logB_m6_vs_flat"] for r in rep_rows], dtype=float)
        arrd = np.array([r["delta_m6_minus_m2"] for r in rep_rows], dtype=float)

        kind_results.append(
            {
                "kind": kind,
                "full_logB_m2_vs_flat": full_m2,
                "full_logB_m6_vs_flat": full_m6,
                "full_delta_m6_minus_m2": full_delta,
                "replications": len(rep_rows),
                "prob_m6_gt_flat": float(np.mean(arr6 > 0.0)),
                "prob_m6_gt_m2": float(np.mean(arrd > 0.0)),
                "median_delta_m6_minus_m2": float(np.median(arrd)),
                "std_delta_m6_minus_m2": float(np.std(arrd)),
                "fit_m6": fit6,
                "resid_m2_rho_theta": rs2["rho_theta_resid"],
                "resid_m6_rho_theta": rs6["rho_theta_resid"],
                "resid_m2_max_abs_bin_mean": rs2["max_abs_bin_mean"],
                "resid_m6_max_abs_bin_mean": rs6["max_abs_bin_mean"],
                "resid_rho_improvement_abs": float(abs(rs2["rho_theta_resid"]) - abs(rs6["rho_theta_resid"])),
                "resid_bin_improvement": float(rs2["max_abs_bin_mean"] - rs6["max_abs_bin_mean"]),
            }
        )

    best = max(
        kind_results,
        key=lambda r: (
            float(r["full_delta_m6_minus_m2"]),
            float(r["prob_m6_gt_m2"]),
            -float(r["std_delta_m6_minus_m2"]),
        ),
    )

    pass_full = best["full_delta_m6_minus_m2"] > 0.0
    pass_rep = best["prob_m6_gt_m2"] >= 0.75 and best["prob_m6_gt_flat"] >= 0.95
    pass_disp = best["std_delta_m6_minus_m2"] <= 0.25
    pass_resid = best["resid_rho_improvement_abs"] > 0.03 and best["resid_bin_improvement"] > 0.03

    if pass_full and pass_rep and pass_disp and pass_resid:
        verdict = "PHYSICAL_DECOHERENCE_EXTENSION_SUPPORTED"
    elif pass_full and pass_rep:
        verdict = "PHYSICAL_DECOHERENCE_EXTENSION_PARTIAL"
    else:
        verdict = "PHYSICAL_DECOHERENCE_EXTENSION_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "locked_protocol_used": {
            "fraction": frac,
            "q_width": q_width,
            "cohort_name": cohort.get("name"),
            "n_pairs": len(rows),
        },
        "kind_results": kind_results,
        "best_kind": best,
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
        "# RAPORT QW-1801: PHYSICAL DECOHERENCE EXTENSION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Locked protocol: frac={frac}, q_width={q_width}, cohort={cohort.get('name')}, n={len(rows)}",
        f"- Best damping kind: {best['kind']}",
        f"- Full delta (M6-M2): {best['full_delta_m6_minus_m2']:.4f}",
        f"- P(M6>M2): {best['prob_m6_gt_m2']:.3f}",
        f"- P(M6>flat): {best['prob_m6_gt_flat']:.3f}",
        f"- Std delta (M6-M2): {best['std_delta_m6_minus_m2']:.3f}",
        f"- Residual improvements: d|rho|={best['resid_rho_improvement_abs']:.4f}, d(bin)={best['resid_bin_improvement']:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_gain_over_m2: {pass_full}",
        f"- replication_gain: {pass_rep}",
        f"- dispersion_control: {pass_disp}",
        f"- residual_flattening: {pass_resid}",
        "",
        "## Kind Table",
    ]
    for r in kind_results:
        lines.append(
            f"- {r['kind']}: full_delta={r['full_delta_m6_minus_m2']:.4f}, P(M6>M2)={r['prob_m6_gt_m2']:.3f}, "
            f"P(M6>flat)={r['prob_m6_gt_flat']:.3f}, std_delta={r['std_delta_m6_minus_m2']:.3f}, "
            f"d|rho|={r['resid_rho_improvement_abs']:.4f}, d(bin)={r['resid_bin_improvement']:.4f}, "
            f"lambda*={r['fit_m6']['lambda']:.4f}"
        )
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1801] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1801] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
