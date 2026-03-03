#!/usr/bin/env python3
"""
QW-1833: GW multi-detector consistency gate on windowed objective space.

Uses QW-1831 features to check whether shared-vs-control signal is
multi-detector consistent under quantile-style diagnostics.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
IN_CSV = ROOT / "gw1831_window_features.csv"
OUT_JSON = ROOT / "report_qw1833_gw_multi_detector_consistency_gate.json"
OUT_MD = ROOT / "RAPORT_QW1833_GW_MULTI_DETECTOR_CONSISTENCY_GATE.md"


def auc_pos_gt_neg(pos: np.ndarray, neg: np.ndarray) -> float:
    """AUC where higher values indicate shared-like score."""
    # rank-based AUC on concatenated arrays
    y = np.concatenate([np.ones(len(pos), dtype=int), np.zeros(len(neg), dtype=int)])
    s = np.concatenate([pos, neg])
    n1 = len(pos)
    n0 = len(neg)
    order = np.argsort(s)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    rs = float(np.sum(ranks[y == 1]))
    return float((rs - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def main() -> None:
    if not IN_CSV.exists():
        raise FileNotFoundError(f"Missing input: {IN_CSV}")

    df = pd.read_csv(IN_CSV)
    req = ["pair", "max_abs_corr", "mean_abs_corr", "corr_at_10ms", "corr_at_0ms"]
    for c in req:
        if c not in df.columns:
            raise RuntimeError(f"Missing column: {c}")

    d = df.dropna(subset=req).copy()

    # Composite shared-likeness score (nonnegative dominant terms + signed lag terms)
    d["score_shared_like"] = (
        0.55 * d["max_abs_corr"].astype(float)
        + 0.25 * d["mean_abs_corr"].astype(float)
        + 0.10 * d["corr_at_0ms"].astype(float)
        + 0.10 * d["corr_at_10ms"].astype(float)
    )

    s_hl = d[d["pair"] == "H1-L1"]["score_shared_like"].to_numpy(dtype=float)
    s_hv = d[d["pair"] == "H1-V1"]["score_shared_like"].to_numpy(dtype=float)
    s_lv = d[d["pair"] == "L1-V1"]["score_shared_like"].to_numpy(dtype=float)
    s_ctrl = np.concatenate([s_hv, s_lv])

    med_hl = float(np.median(s_hl))
    med_hv = float(np.median(s_hv))
    med_lv = float(np.median(s_lv))
    med_ctrl = float(np.median(s_ctrl))

    # separation and consistency metrics
    sep_vs_ctrl = float(med_hl - med_ctrl)
    ctrl_gap = float(abs(med_hv - med_lv))
    ctrl_std = float(np.std([med_hv, med_lv]))

    q90_ctrl = float(np.quantile(s_ctrl, 0.90))
    p_shared_above_q90 = float(np.mean(s_hl > q90_ctrl))
    p_ctrl_above_q90 = float(np.mean(s_ctrl > q90_ctrl))
    adv_q90 = float(p_shared_above_q90 - p_ctrl_above_q90)

    auc_hl_vs_hv = auc_pos_gt_neg(s_hl, s_hv)
    auc_hl_vs_lv = auc_pos_gt_neg(s_hl, s_lv)
    auc_hl_vs_ctrl = auc_pos_gt_neg(s_hl, s_ctrl)

    summary = {
        "n_h1l1": int(len(s_hl)),
        "n_h1v1": int(len(s_hv)),
        "n_l1v1": int(len(s_lv)),
        "median_score_h1l1": med_hl,
        "median_score_h1v1": med_hv,
        "median_score_l1v1": med_lv,
        "median_score_ctrl": med_ctrl,
        "sep_median_h1l1_minus_ctrl": sep_vs_ctrl,
        "control_median_gap": ctrl_gap,
        "control_median_std": ctrl_std,
        "p_shared_above_ctrl_q90": p_shared_above_q90,
        "p_ctrl_above_ctrl_q90": p_ctrl_above_q90,
        "adv_shared_minus_ctrl_q90": adv_q90,
        "auc_h1l1_vs_h1v1": float(auc_hl_vs_hv),
        "auc_h1l1_vs_l1v1": float(auc_hl_vs_lv),
        "auc_h1l1_vs_ctrl_mix": float(auc_hl_vs_ctrl),
    }

    pass_sep = summary["sep_median_h1l1_minus_ctrl"] >= 0.002
    pass_adv = summary["adv_shared_minus_ctrl_q90"] >= 0.30
    pass_auc = summary["auc_h1l1_vs_ctrl_mix"] >= 0.75
    pass_ctrl_cons = summary["control_median_gap"] <= 0.0025

    if pass_sep and pass_adv and pass_auc and pass_ctrl_cons:
        verdict = "GW_MULTI_DETECTOR_CONSISTENT"
    elif pass_sep and pass_adv and pass_auc:
        verdict = "GW_MULTI_DETECTOR_PARTIAL_CONTROL_MISMATCH"
    else:
        verdict = "GW_MULTI_DETECTOR_INCONSISTENT"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "pass_flags": {
            "shared_vs_ctrl_separation": bool(pass_sep),
            "shared_prevalence_advantage": bool(pass_adv),
            "auc_support": bool(pass_auc),
            "control_pair_consistency": bool(pass_ctrl_cons),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1833: GW MULTI-DETECTOR CONSISTENCY GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- n(H1-L1/H1-V1/L1-V1): {summary['n_h1l1']}/{summary['n_h1v1']}/{summary['n_l1v1']}",
        f"- Median separation H1-L1 - CTRL: {summary['sep_median_h1l1_minus_ctrl']:.6f}",
        f"- Control median gap (H1-V1 vs L1-V1): {summary['control_median_gap']:.6f}",
        f"- Advantage above CTRL q90: {summary['adv_shared_minus_ctrl_q90']:.3f}",
        f"- AUC H1-L1 vs CTRL mix: {summary['auc_h1l1_vs_ctrl_mix']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- shared_vs_ctrl_separation: {pass_sep}",
        f"- shared_prevalence_advantage: {pass_adv}",
        f"- auc_support: {pass_auc}",
        f"- control_pair_consistency: {pass_ctrl_cons}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1833] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1833] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
