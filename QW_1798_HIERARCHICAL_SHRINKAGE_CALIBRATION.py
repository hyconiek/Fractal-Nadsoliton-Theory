#!/usr/bin/env python3
"""
QW-1798: Hierarchical shrinkage calibration.

Follow-up to QW-1797:
- QW-1797 showed strong gain but high replication variance.
- Here we calibrate shrinkage strength to improve stability while keeping gain.
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
OUT_JSON = ROOT / "report_qw1798_hierarchical_shrinkage_calibration.json"
OUT_MD = ROOT / "RAPORT_QW1798_HIERARCHICAL_SHRINKAGE_CALIBRATION.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def main() -> None:
    helper = load_module(ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py", "qw1787_helper")
    core = load_module(ROOT / "QW_1797_HIERARCHICAL_MULTIMODE_SHRINKAGE.py", "qw1797_core")

    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))
    d1796 = json.loads((ROOT / "report_qw1796_physical_multimode_extension.json").read_text(encoding="utf-8"))

    q_center = float(d1773["projection"]["p"])
    q_width = float(d1793["operational_protocol"]["q_width"])
    frac = float(d1793["operational_protocol"]["fraction"])
    cohort = d1793["operational_protocol"]["cohort"]
    n_match_min = int(cohort["n_match_min"])
    stability_max = float(cohort["stability_max"])
    ell_list = list(d1796["best_family"]["ell_list"])
    family_name = str(d1796["best_family"]["family"])

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

    raw_prior = core.posterior_hyper_from_full(
        helper, theta_all, H_all, q_center=q_center, q_width=q_width, ell_list=ell_list, n_mc=60000, seed=17981
    )
    raw_b_std = float(raw_prior["posterior_b_std_raw"])
    raw_a_std = float(raw_prior["posterior_alpha_std_raw"])

    z0 = core.evidence_flat(H_all, n_mc=8000, seed=17982)
    z2 = core.evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=12000, seed=17983)
    full_m2 = float(z2 - z0)

    shrink_factors = [0.20, 0.35, 0.50, 0.70, 1.00]
    sweep = []
    for j, sf in enumerate(shrink_factors):
        prior = {
            "b_mean": raw_prior["b_mean"],
            "alpha_mean": raw_prior["alpha_mean"],
            "b_std": max(0.03, sf * raw_b_std),
            "alpha_std": max(0.05, sf * raw_a_std),
        }

        z5 = core.evidence_m5(
            helper,
            theta_all,
            H_all,
            q_center=q_center,
            q_width=q_width,
            ell_list=ell_list,
            n_mc=15000,
            seed=18000 + 100 * j + 1,
            b_prior_mean=prior["b_mean"],
            b_prior_std=prior["b_std"],
            alpha_prior_mean=prior["alpha_mean"],
            alpha_prior_std=prior["alpha_std"],
        )
        full_m5 = float(z5 - z0)
        full_delta = float(full_m5 - full_m2)

        rng = np.random.default_rng(18000 + 100 * j + 2)
        rep_rows = []
        for i in range(12):
            idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
            if len(idx) < 80:
                continue
            th = theta_all[idx]
            hh = H_all[idx]
            b0 = core.evidence_flat(hh, n_mc=3800, seed=18000 + 100 * j + 20 * i + 3)
            b2 = core.evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=18000 + 100 * j + 20 * i + 4)
            b5 = core.evidence_m5(
                helper,
                th,
                hh,
                q_center=q_center,
                q_width=q_width,
                ell_list=ell_list,
                n_mc=8200,
                seed=18000 + 100 * j + 20 * i + 5,
                b_prior_mean=prior["b_mean"],
                b_prior_std=prior["b_std"],
                alpha_prior_mean=prior["alpha_mean"],
                alpha_prior_std=prior["alpha_std"],
            )
            l2 = float(b2 - b0)
            l5 = float(b5 - b0)
            rep_rows.append({"delta_hier_vs_m2": float(l5 - l2), "logB_hier_vs_flat": l5})

        arr_d = np.array([r["delta_hier_vs_m2"] for r in rep_rows], dtype=float)
        arr_h = np.array([r["logB_hier_vs_flat"] for r in rep_rows], dtype=float)
        p_gain = float(np.mean(arr_d > 0.0))
        p_pos = float(np.mean(arr_h > 0.0))
        std_d = float(np.std(arr_d))
        med_d = float(np.median(arr_d))

        score = (
            0.35 * clip01(full_delta / 1.5)
            + 0.25 * p_gain
            + 0.20 * p_pos
            + 0.20 * clip01(1.0 - std_d / 0.40)
        )

        pass_basic = full_delta > 0.0 and p_gain >= 0.80 and p_pos >= 0.95 and std_d <= 0.35

        sweep.append(
            {
                "shrink_factor": sf,
                "b_std": prior["b_std"],
                "alpha_std": prior["alpha_std"],
                "full_logB_m2_vs_flat": full_m2,
                "full_logB_hier_vs_flat": full_m5,
                "full_delta_hier_vs_m2": full_delta,
                "replications": len(rep_rows),
                "prob_hier_gt_m2": p_gain,
                "prob_hier_gt_flat": p_pos,
                "median_delta_hier_vs_m2": med_d,
                "std_delta_hier_vs_m2": std_d,
                "selection_score": float(score),
                "pass_basic": bool(pass_basic),
            }
        )

    valid = [r for r in sweep if r["pass_basic"]]
    if len(valid) > 0:
        best = max(valid, key=lambda r: r["selection_score"])
        recommendation_strength = "STRONG"
    else:
        best = max(sweep, key=lambda r: r["selection_score"])
        recommendation_strength = "CONDITIONAL"

    # Residual diagnostic for best shrink.
    best_prior = {
        "b_mean": raw_prior["b_mean"],
        "alpha_mean": raw_prior["alpha_mean"],
        "b_std": float(best["b_std"]),
        "alpha_std": float(best["alpha_std"]),
    }
    fit2 = core.fit_best_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=65000, seed=18901)
    fit5 = core.fit_best_m5_hier(
        helper, theta_all, H_all, q_center=q_center, q_width=q_width, ell_list=ell_list, prior=best_prior, n_mc=85000, seed=18902
    )
    hd0 = np.clip(helper.hellings_downs(theta_all), 1e-9, None)
    pred2 = fit2["A"] * (hd0 ** fit2["q"]) + fit2["C"]
    T5 = core.template_family(theta_all, ell_list, fit5["alpha"])
    pred5 = fit5["A"] * (hd0 ** fit5["q"]) + fit5["B"] * T5 + fit5["C"]
    rs2 = core.residual_stats(H_all - pred2, theta_all)
    rs5 = core.residual_stats(H_all - pred5, theta_all)
    residual_improvement = {
        "d_abs_rho": float(abs(rs2["rho_theta_resid"]) - abs(rs5["rho_theta_resid"])),
        "d_bin_flatness": float(rs2["max_abs_bin_mean"] - rs5["max_abs_bin_mean"]),
    }

    if recommendation_strength == "STRONG":
        verdict = "HIERARCHICAL_SHRINKAGE_CALIBRATED"
    else:
        verdict = "HIERARCHICAL_SHRINKAGE_PARTIAL"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "family": family_name,
        "ell_list": ell_list,
        "n_pairs": len(rows),
        "fraction": frac,
        "q_width": q_width,
        "raw_posterior_prior": raw_prior,
        "sweep": sweep,
        "best": best,
        "recommendation_strength": recommendation_strength,
        "residual_improvement_best": residual_improvement,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1798: HIERARCHICAL SHRINKAGE CALIBRATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Family: {family_name} {ell_list}",
        f"- Best shrink factor: {best['shrink_factor']:.2f} ({recommendation_strength})",
        f"- Full delta hier-M2 (best): {best['full_delta_hier_vs_m2']:.4f}",
        f"- P(hier>M2) (best): {best['prob_hier_gt_m2']:.3f}",
        f"- P(hier>flat) (best): {best['prob_hier_gt_flat']:.3f}",
        f"- Std delta hier-M2 (best): {best['std_delta_hier_vs_m2']:.3f}",
        f"- Residual improvements (best): d|rho|={residual_improvement['d_abs_rho']:.4f}, d(bin)={residual_improvement['d_bin_flatness']:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Sweep",
    ]
    for r in sweep:
        lines.append(
            f"- sf={r['shrink_factor']:.2f}: full_delta={r['full_delta_hier_vs_m2']:.4f}, "
            f"P(hier>M2)={r['prob_hier_gt_m2']:.3f}, P(hier>flat)={r['prob_hier_gt_flat']:.3f}, "
            f"std_delta={r['std_delta_hier_vs_m2']:.3f}, score={r['selection_score']:.3f}, pass_basic={r['pass_basic']}"
        )
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1798] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1798] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
