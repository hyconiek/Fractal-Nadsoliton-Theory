#!/usr/bin/env python3
"""
QW-1804: Variance decomposition for additive-hierarchical branch.

Question:
- Is high replication dispersion mainly due to Monte Carlo evidence estimator noise
  or due to true split-to-split data instability?
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1804_variance_decomposition_mc.json"
OUT_MD = ROOT / "RAPORT_QW1804_VARIANCE_DECOMPOSITION_MC.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def stratified_split(theta: np.ndarray, rng: np.random.Generator, test_frac: float = 0.25, n_bins: int = 8) -> Tuple[np.ndarray, np.ndarray]:
    bins = np.linspace(0.0, 180.0, n_bins + 1)
    idx_all = np.arange(len(theta))
    test = []
    for i in range(n_bins):
        m = (theta >= bins[i]) & (theta < bins[i + 1] if i < n_bins - 1 else theta <= bins[i + 1])
        idx = idx_all[m]
        if len(idx) == 0:
            continue
        k = max(1, int(round(test_frac * len(idx))))
        k = min(k, len(idx) - 1) if len(idx) > 1 else 1
        test.append(rng.choice(idx, size=k, replace=False))
    if len(test) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    test_idx = np.sort(np.unique(np.concatenate(test)))
    train_mask = np.ones(len(theta), dtype=bool)
    train_mask[test_idx] = False
    train_idx = np.where(train_mask)[0]
    return train_idx, test_idx


def main() -> None:
    helper = load_module(ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py", "qw1787_helper")
    addcore = load_module(ROOT / "QW_1803_ADDITIVE_HIERARCHICAL_STABILIZATION.py", "qw1803_add")

    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))
    d1803 = json.loads((ROOT / "report_qw1803_additive_hierarchical_stabilization.json").read_text(encoding="utf-8"))

    q_center = float(d1773["projection"]["p"])
    q_width = float(d1793["operational_protocol"]["q_width"])
    cohort = d1793["operational_protocol"]["cohort"]
    n_match_min = int(cohort["n_match_min"])
    stability_max = float(cohort["stability_max"])

    raw = d1803["raw_posterior"]
    best = d1803["best"]
    prior = {
        "b_mean": float(raw["b_mean"]),
        "lam_mean": float(raw["lam_mean"]),
        "b_std": float(best["b_std"]),
        "lam_std": float(best["lam_std"]),
    }

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

    rng = np.random.default_rng(18040)
    split_count = 6
    mc_repeats = 6
    split_records = []
    split_means = []
    split_within_vars = []

    for s in range(split_count):
        tr_idx, te_idx = stratified_split(theta_all, rng=rng, test_frac=0.25, n_bins=8)
        if len(tr_idx) < 50:
            continue
        th = theta_all[tr_idx]
        hh = H_all[tr_idx]

        deltas = []
        for r in range(mc_repeats):
            seed_base = 19000 + 1000 * s + 50 * r
            b0 = addcore.evidence_flat(hh, n_mc=3800, seed=seed_base + 1)
            b2 = addcore.evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=seed_base + 2)
            b6 = addcore.evidence_m6_add(
                helper,
                th,
                hh,
                q_center=q_center,
                q_width=q_width,
                n_mc=8200,
                seed=seed_base + 3,
                b_mean=prior["b_mean"],
                b_std=prior["b_std"],
                lam_mean=prior["lam_mean"],
                lam_std=prior["lam_std"],
            )
            deltas.append(float((b6 - b0) - (b2 - b0)))

        arr = np.array(deltas, dtype=float)
        mean_delta = float(np.mean(arr))
        var_within = float(np.var(arr, ddof=1)) if len(arr) > 1 else 0.0
        split_records.append(
            {
                "split": s,
                "n_train": int(len(tr_idx)),
                "mc_repeats": mc_repeats,
                "delta_samples": deltas,
                "mean_delta": mean_delta,
                "std_delta_mc": float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
            }
        )
        split_means.append(mean_delta)
        split_within_vars.append(var_within)

    means = np.array(split_means, dtype=float)
    within = float(np.mean(split_within_vars)) if len(split_within_vars) > 0 else 0.0
    between = float(np.var(means, ddof=1)) if len(means) > 1 else 0.0
    total = within + between
    within_ratio = float(within / total) if total > 0 else 0.0
    between_ratio = float(between / total) if total > 0 else 0.0

    summary = {
        "n_pairs": len(rows),
        "split_count": len(split_records),
        "mc_repeats_per_split": mc_repeats,
        "mean_delta_overall": float(np.mean(means)) if len(means) > 0 else 0.0,
        "median_delta_overall": float(np.median(means)) if len(means) > 0 else 0.0,
        "std_of_split_means": float(np.std(means, ddof=1)) if len(means) > 1 else 0.0,
        "within_variance_mc": within,
        "between_variance_splits": between,
        "within_ratio": within_ratio,
        "between_ratio": between_ratio,
    }

    pass_positive = summary["mean_delta_overall"] > 0.0
    pass_mc_dominant = within_ratio >= 0.55
    pass_between_small = summary["std_of_split_means"] <= 0.50

    if pass_positive and pass_mc_dominant and pass_between_small:
        verdict = "VARIANCE_MAINLY_ESTIMATOR_DRIVEN"
    elif pass_positive and pass_between_small:
        verdict = "VARIANCE_MIXED_BUT_MANAGEABLE"
    else:
        verdict = "VARIANCE_MAINLY_DATA_INSTABILITY"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "prior_used": prior,
        "summary": summary,
        "pass_flags": {
            "positive_mean_delta": bool(pass_positive),
            "mc_dominant_component": bool(pass_mc_dominant),
            "between_split_std_small": bool(pass_between_small),
        },
        "verdict": verdict,
        "split_records": split_records,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1804: VARIANCE DECOMPOSITION MC",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Splits: {summary['split_count']}, MC repeats per split: {mc_repeats}",
        f"- Mean delta overall: {summary['mean_delta_overall']:.4f}",
        f"- Std of split means: {summary['std_of_split_means']:.4f}",
        f"- Within variance (MC): {within:.6f}",
        f"- Between variance (splits): {between:.6f}",
        f"- Variance ratio within/between: {within_ratio:.3f}/{between_ratio:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- positive_mean_delta: {pass_positive}",
        f"- mc_dominant_component: {pass_mc_dominant}",
        f"- between_split_std_small: {pass_between_small}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1804] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1804] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
