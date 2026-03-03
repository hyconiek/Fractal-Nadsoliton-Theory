#!/usr/bin/env python3
"""
QW-1790: Residual angular augmentation test.

Model extension:
- M2: reparam HD^q + C
- M3: reparam HD^q + B*P2(cos(theta)) + C
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
OUT_JSON = ROOT / "report_qw1790_residual_angular_augmentation.json"
OUT_MD = ROOT / "RAPORT_QW1790_RESIDUAL_ANGULAR_AUGMENTATION.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def evidence_reparam_width(helper, theta: np.ndarray, H: np.ndarray, sigma: float, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array(
        [helper.loglike(theta, H, sigma, a * (hd0 ** qq) + c) for a, c, qq in zip(A, C, q)],
        dtype=float,
    )
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_augmented(helper, theta: np.ndarray, H: np.ndarray, sigma: float, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B = rng.uniform(-1.0, 1.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    th = np.radians(theta)
    p2 = 0.5 * (3.0 * np.cos(th) ** 2 - 1.0)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array(
        [helper.loglike(theta, H, sigma, a * (hd0 ** qq) + b * p2 + c) for a, b, c, qq in zip(A, B, C, q)],
        dtype=float,
    )
    return float(logsumexp(ll) - np.log(len(ll)))


def main() -> None:
    helper = load_helper()
    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    q_center = float(d1773["projection"]["p"])
    q_width = 0.20
    frac = 0.95
    n_rep = 16

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
        if len(x) >= 120 and stab <= 0.65:
            rows.append({"theta_deg": float(sep), "hxy": float(hxy)})

    if len(rows) < 85:
        raise RuntimeError("Base operational cohort too small.")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    sigma_all = max(float(np.std(H_all)), 1e-6)

    lz0 = helper.evidence_flat(theta_all, H_all, sigma_all, n_mc=10000, seed=9001)
    lz1 = helper.evidence_legacy(theta_all, H_all, sigma_all, n_mc=13000, seed=9002)
    lz2 = evidence_reparam_width(helper, theta_all, H_all, sigma_all, q_center=q_center, q_width=q_width, n_mc=17000, seed=9003)
    lz3 = evidence_augmented(helper, theta_all, H_all, sigma_all, q_center=q_center, q_width=q_width, n_mc=23000, seed=9004)

    full_b1 = float(lz1 - lz0)
    full_b2 = float(lz2 - lz0)
    full_b3 = float(lz3 - lz0)
    full_delta_23 = float(full_b3 - full_b2)

    rng = np.random.default_rng(1790)
    rep_rows = []
    for i in range(n_rep):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 80:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        sg = max(float(np.std(hh)), 1e-6)
        b0 = helper.evidence_flat(th, hh, sg, n_mc=6000, seed=9100 + 10 * i + 1)
        b2 = evidence_reparam_width(helper, th, hh, sg, q_center=q_center, q_width=q_width, n_mc=9000, seed=9100 + 10 * i + 2)
        b3 = evidence_augmented(helper, th, hh, sg, q_center=q_center, q_width=q_width, n_mc=12000, seed=9100 + 10 * i + 3)
        rep_rows.append(
            {
                "rep": i,
                "n_pairs": int(len(idx)),
                "logB_m2_reparam_vs_flat": float(b2 - b0),
                "logB_m3_aug_vs_flat": float(b3 - b0),
                "delta_logB_m3_minus_m2": float((b3 - b0) - (b2 - b0)),
            }
        )

    arr_m2 = np.array([r["logB_m2_reparam_vs_flat"] for r in rep_rows], dtype=float)
    arr_m3 = np.array([r["logB_m3_aug_vs_flat"] for r in rep_rows], dtype=float)
    arr_d = np.array([r["delta_logB_m3_minus_m2"] for r in rep_rows], dtype=float)

    summary = {
        "n_pairs_operational": len(rows),
        "q_width": q_width,
        "fraction": frac,
        "replications": len(rep_rows),
        "full_logB_legacy_vs_flat": full_b1,
        "full_logB_m2_reparam_vs_flat": full_b2,
        "full_logB_m3_aug_vs_flat": full_b3,
        "full_delta_logB_m3_minus_m2": full_delta_23,
        "prob_m3_gt_m2": float(np.mean(arr_d > 0.0)),
        "prob_m3_gt_flat": float(np.mean(arr_m3 > 0.0)),
        "median_delta_logB_m3_minus_m2": float(np.median(arr_d)),
        "std_delta_logB_m3_minus_m2": float(np.std(arr_d)),
    }

    pass_full_improvement = summary["full_delta_logB_m3_minus_m2"] > 0.0
    pass_rep_improvement = summary["prob_m3_gt_m2"] >= 0.80
    pass_rep_positive = summary["prob_m3_gt_flat"] >= 0.95
    pass_dispersion = summary["std_delta_logB_m3_minus_m2"] <= 0.18

    if pass_full_improvement and pass_rep_improvement and pass_rep_positive and pass_dispersion:
        verdict = "ANGULAR_AUGMENTATION_SUPPORTED"
    elif pass_full_improvement and pass_rep_improvement:
        verdict = "ANGULAR_AUGMENTATION_PARTIAL"
    else:
        verdict = "ANGULAR_AUGMENTATION_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "pass_flags": {
            "full_improvement": bool(pass_full_improvement),
            "replication_improvement": bool(pass_rep_improvement),
            "replication_positive_m3": bool(pass_rep_positive),
            "dispersion_control": bool(pass_dispersion),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1790: RESIDUAL ANGULAR AUGMENTATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Operational pairs: {len(rows)}",
        f"- Full logB M2 vs flat: {full_b2:.4f}",
        f"- Full logB M3 vs flat: {full_b3:.4f}",
        f"- Full delta (M3-M2): {full_delta_23:.4f}",
        f"- P(M3>M2): {summary['prob_m3_gt_m2']:.3f}",
        f"- P(M3>flat): {summary['prob_m3_gt_flat']:.3f}",
        f"- Std delta (M3-M2): {summary['std_delta_logB_m3_minus_m2']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_improvement: {pass_full_improvement}",
        f"- replication_improvement: {pass_rep_improvement}",
        f"- replication_positive_m3: {pass_rep_positive}",
        f"- dispersion_control: {pass_dispersion}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1790] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1790] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
