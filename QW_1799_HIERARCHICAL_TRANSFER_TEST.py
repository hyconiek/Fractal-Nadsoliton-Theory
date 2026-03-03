#!/usr/bin/env python3
"""
QW-1799: Hierarchical prior transfer test (train->test).

Goal:
- test whether hierarchical multimode gain survives on disjoint holdout sets,
- detect overfitting / leakage-like behavior in the strong-gain regime from QW-1797/1798.
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
OUT_JSON = ROOT / "report_qw1799_hierarchical_transfer_test.json"
OUT_MD = ROOT / "RAPORT_QW1799_HIERARCHICAL_TRANSFER_TEST.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def stratified_split(theta: np.ndarray, rng: np.random.Generator, test_frac: float = 0.45, n_bins: int = 8) -> tuple[np.ndarray, np.ndarray]:
    bins = np.linspace(0.0, 180.0, n_bins + 1)
    idx_all = np.arange(len(theta))
    test_idx = []
    for i in range(n_bins):
        m = (theta >= bins[i]) & (theta < bins[i + 1] if i < n_bins - 1 else theta <= bins[i + 1])
        idx = idx_all[m]
        if len(idx) == 0:
            continue
        k = max(1, int(round(test_frac * len(idx))))
        k = min(k, len(idx) - 1) if len(idx) > 1 else 1
        test_idx.append(rng.choice(idx, size=k, replace=False))
    if len(test_idx) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    test_idx = np.sort(np.unique(np.concatenate(test_idx)))
    train_mask = np.ones(len(theta), dtype=bool)
    train_mask[test_idx] = False
    train_idx = np.where(train_mask)[0]
    return train_idx, test_idx


def main() -> None:
    helper = load_module(ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py", "qw1787_helper")
    core = load_module(ROOT / "QW_1797_HIERARCHICAL_MULTIMODE_SHRINKAGE.py", "qw1797_core")

    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))
    d1796 = json.loads((ROOT / "report_qw1796_physical_multimode_extension.json").read_text(encoding="utf-8"))
    d1798 = json.loads((ROOT / "report_qw1798_hierarchical_shrinkage_calibration.json").read_text(encoding="utf-8"))

    q_center = float(d1773["projection"]["p"])
    q_width = float(d1793["operational_protocol"]["q_width"])
    cohort = d1793["operational_protocol"]["cohort"]
    n_match_min = int(cohort["n_match_min"])
    stability_max = float(cohort["stability_max"])
    ell_list = list(d1796["best_family"]["ell_list"])
    family_name = str(d1796["best_family"]["family"])
    best_shrink = float(d1798["best"]["shrink_factor"])

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

    # Global raw posterior moments (for initial scale only).
    raw = core.posterior_hyper_from_full(
        helper, theta_all, H_all, q_center=q_center, q_width=q_width, ell_list=ell_list, n_mc=50000, seed=17991
    )
    base_b_std = max(0.03, best_shrink * float(raw["posterior_b_std_raw"]))
    base_a_std = max(0.05, best_shrink * float(raw["posterior_alpha_std_raw"]))

    rng = np.random.default_rng(17999)
    split_rows = []
    for rep in range(14):
        tr_idx, te_idx = stratified_split(theta_all, rng=rng, test_frac=0.45, n_bins=8)
        if len(tr_idx) < 35 or len(te_idx) < 35:
            continue

        th_tr = theta_all[tr_idx]
        hh_tr = H_all[tr_idx]
        th_te = theta_all[te_idx]
        hh_te = H_all[te_idx]

        # Learn hierarchical means on TRAIN.
        post_tr = core.posterior_hyper_from_full(
            helper, th_tr, hh_tr, q_center=q_center, q_width=q_width, ell_list=ell_list, n_mc=35000, seed=18000 + 100 * rep + 1
        )
        prior = {
            "b_mean": post_tr["b_mean"],
            "alpha_mean": post_tr["alpha_mean"],
            "b_std": base_b_std,
            "alpha_std": base_a_std,
        }

        # TRAIN evidences
        z0_tr = core.evidence_flat(hh_tr, n_mc=3200, seed=18000 + 100 * rep + 2)
        z2_tr = core.evidence_m2(helper, th_tr, hh_tr, q_center=q_center, q_width=q_width, n_mc=5200, seed=18000 + 100 * rep + 3)
        z5h_tr = core.evidence_m5(
            helper, th_tr, hh_tr, q_center=q_center, q_width=q_width, ell_list=ell_list, n_mc=7000, seed=18000 + 100 * rep + 4,
            b_prior_mean=prior["b_mean"], b_prior_std=prior["b_std"], alpha_prior_mean=prior["alpha_mean"], alpha_prior_std=prior["alpha_std"]
        )
        delta_tr = float((z5h_tr - z0_tr) - (z2_tr - z0_tr))

        # TEST evidences (same learned prior)
        z0_te = core.evidence_flat(hh_te, n_mc=3200, seed=18000 + 100 * rep + 5)
        z2_te = core.evidence_m2(helper, th_te, hh_te, q_center=q_center, q_width=q_width, n_mc=5200, seed=18000 + 100 * rep + 6)
        z5h_te = core.evidence_m5(
            helper, th_te, hh_te, q_center=q_center, q_width=q_width, ell_list=ell_list, n_mc=7000, seed=18000 + 100 * rep + 7,
            b_prior_mean=prior["b_mean"], b_prior_std=prior["b_std"], alpha_prior_mean=prior["alpha_mean"], alpha_prior_std=prior["alpha_std"]
        )
        l2_te = float(z2_te - z0_te)
        l5_te = float(z5h_te - z0_te)
        delta_te = float(l5_te - l2_te)

        split_rows.append(
            {
                "rep": rep,
                "n_train": int(len(tr_idx)),
                "n_test": int(len(te_idx)),
                "train_delta_hier_vs_m2": delta_tr,
                "test_logB_m2_vs_flat": l2_te,
                "test_logB_hier_vs_flat": l5_te,
                "test_delta_hier_vs_m2": delta_te,
                "generalization_gap_delta": float(delta_tr - delta_te),
            }
        )

    arr_te = np.array([r["test_delta_hier_vs_m2"] for r in split_rows], dtype=float)
    arr_gap = np.array([r["generalization_gap_delta"] for r in split_rows], dtype=float)
    arr_hf = np.array([r["test_logB_hier_vs_flat"] for r in split_rows], dtype=float)

    summary = {
        "family": family_name,
        "ell_list": ell_list,
        "n_pairs_total": len(rows),
        "shrink_factor": best_shrink,
        "prior_std": {"b_std": base_b_std, "alpha_std": base_a_std},
        "splits": len(split_rows),
        "test_prob_hier_gt_m2": float(np.mean(arr_te > 0.0)),
        "test_prob_hier_gt_flat": float(np.mean(arr_hf > 0.0)),
        "test_median_delta_hier_vs_m2": float(np.median(arr_te)),
        "test_std_delta_hier_vs_m2": float(np.std(arr_te)),
        "gap_median_train_minus_test": float(np.median(arr_gap)),
        "gap_std_train_minus_test": float(np.std(arr_gap)),
    }

    pass_transfer_gain = summary["test_prob_hier_gt_m2"] >= 0.75 and summary["test_median_delta_hier_vs_m2"] > 0.0
    pass_transfer_stability = summary["test_std_delta_hier_vs_m2"] <= 0.35
    pass_gap = summary["gap_median_train_minus_test"] <= 0.90
    pass_flat = summary["test_prob_hier_gt_flat"] >= 0.90

    if pass_transfer_gain and pass_transfer_stability and pass_gap and pass_flat:
        verdict = "HIERARCHICAL_TRANSFER_SUPPORTED"
    elif pass_transfer_gain and pass_flat:
        verdict = "HIERARCHICAL_TRANSFER_PARTIAL"
    else:
        verdict = "HIERARCHICAL_TRANSFER_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "pass_flags": {
            "transfer_gain": bool(pass_transfer_gain),
            "transfer_stability": bool(pass_transfer_stability),
            "generalization_gap_control": bool(pass_gap),
            "test_positive_vs_flat": bool(pass_flat),
        },
        "verdict": verdict,
        "split_results": split_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1799: HIERARCHICAL TRANSFER TEST",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Family: {family_name} {ell_list}",
        f"- Shrink factor: {best_shrink:.2f}",
        f"- Splits: {summary['splits']}",
        f"- Test P(hier>M2): {summary['test_prob_hier_gt_m2']:.3f}",
        f"- Test P(hier>flat): {summary['test_prob_hier_gt_flat']:.3f}",
        f"- Test median delta hier-M2: {summary['test_median_delta_hier_vs_m2']:.4f}",
        f"- Test std delta hier-M2: {summary['test_std_delta_hier_vs_m2']:.4f}",
        f"- Gap median (train-test): {summary['gap_median_train_minus_test']:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- transfer_gain: {pass_transfer_gain}",
        f"- transfer_stability: {pass_transfer_stability}",
        f"- generalization_gap_control: {pass_gap}",
        f"- test_positive_vs_flat: {pass_flat}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1799] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1799] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
