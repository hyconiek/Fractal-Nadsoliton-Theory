#!/usr/bin/env python3
"""
QW-1836: GW control-calibrated objective (phase-2 structural reparam).

Idea:
- keep shared score for H1-L1,
- calibrate nuisance offsets for control pairs (H1-V1, L1-V1)
  using training folds only,
- evaluate on blocked test folds.

Goal:
- reduce control-pair mismatch without sacrificing shared-vs-control separation.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
IN_CSV = ROOT / "gw1831_window_features.csv"
OUT_JSON = ROOT / "report_qw1836_gw_control_calibrated_objective.json"
OUT_MD = ROOT / "RAPORT_QW1836_GW_CONTROL_CALIBRATED_OBJECTIVE.md"


def auc_pos_gt_neg(pos: np.ndarray, neg: np.ndarray) -> float:
    y = np.concatenate([np.ones(len(pos), dtype=int), np.zeros(len(neg), dtype=int)])
    s = np.concatenate([pos, neg])
    n1 = len(pos)
    n0 = len(neg)
    order = np.argsort(s)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    rs = float(np.sum(ranks[y == 1]))
    return float((rs - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def score_raw(df: pd.DataFrame) -> np.ndarray:
    return (
        0.55 * df["max_abs_corr"].to_numpy(dtype=float)
        + 0.25 * df["mean_abs_corr"].to_numpy(dtype=float)
        + 0.10 * df["corr_at_0ms"].to_numpy(dtype=float)
        + 0.10 * df["corr_at_10ms"].to_numpy(dtype=float)
    )


def metrics_from_scores(pair: np.ndarray, s: np.ndarray) -> Dict[str, float]:
    sh = s[pair == "H1-L1"]
    hv = s[pair == "H1-V1"]
    lv = s[pair == "L1-V1"]
    ctrl = np.concatenate([hv, lv])

    auc = auc_pos_gt_neg(sh, ctrl)
    q90 = float(np.quantile(ctrl, 0.90))
    p_shared = float(np.mean(sh > q90))
    p_ctrl = float(np.mean(ctrl > q90))
    adv = float(p_shared - p_ctrl)

    return {
        "auc_h1l1_vs_ctrl": float(auc),
        "adv_shared_minus_ctrl_q90": adv,
        "control_gap": float(abs(np.median(hv) - np.median(lv))),
        "p_shared_above_q90": p_shared,
        "p_ctrl_above_q90": p_ctrl,
        "median_h1l1": float(np.median(sh)),
        "median_h1v1": float(np.median(hv)),
        "median_l1v1": float(np.median(lv)),
    }


def main() -> None:
    if not IN_CSV.exists():
        raise FileNotFoundError(f"Missing input CSV: {IN_CSV}")

    df = pd.read_csv(IN_CSV)
    req = ["pair", "window_idx", "max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]
    for c in req:
        if c not in df.columns:
            raise RuntimeError(f"Missing column: {c}")

    d = df.dropna(subset=req).copy()
    d["score_raw"] = score_raw(d)

    fold_rows: List[Dict[str, float]] = []
    n_folds = 5

    for fold in range(n_folds):
        tr = d[(d["window_idx"].astype(int) % n_folds) != fold].copy()
        te = d[(d["window_idx"].astype(int) % n_folds) == fold].copy()

        # Train-time nuisance offsets for control pairs only
        med_hv = float(np.median(tr.loc[tr["pair"] == "H1-V1", "score_raw"]))
        med_lv = float(np.median(tr.loc[tr["pair"] == "L1-V1", "score_raw"]))
        med_ctrl = float(np.median(tr.loc[tr["pair"] != "H1-L1", "score_raw"]))

        off_hv = med_hv - med_ctrl
        off_lv = med_lv - med_ctrl

        te["score_cal"] = te["score_raw"].to_numpy(dtype=float)
        te.loc[te["pair"] == "H1-V1", "score_cal"] = te.loc[te["pair"] == "H1-V1", "score_cal"] - off_hv
        te.loc[te["pair"] == "L1-V1", "score_cal"] = te.loc[te["pair"] == "L1-V1", "score_cal"] - off_lv
        # H1-L1 left unchanged intentionally

        raw_m = metrics_from_scores(te["pair"].to_numpy(), te["score_raw"].to_numpy(dtype=float))
        cal_m = metrics_from_scores(te["pair"].to_numpy(), te["score_cal"].to_numpy(dtype=float))

        fold_rows.append(
            {
                "fold": fold,
                "n_train": int(len(tr)),
                "n_test": int(len(te)),
                "offset_h1v1": off_hv,
                "offset_l1v1": off_lv,
                "raw_auc": raw_m["auc_h1l1_vs_ctrl"],
                "raw_adv": raw_m["adv_shared_minus_ctrl_q90"],
                "raw_control_gap": raw_m["control_gap"],
                "cal_auc": cal_m["auc_h1l1_vs_ctrl"],
                "cal_adv": cal_m["adv_shared_minus_ctrl_q90"],
                "cal_control_gap": cal_m["control_gap"],
            }
        )

    arr_raw_auc = np.array([r["raw_auc"] for r in fold_rows], dtype=float)
    arr_raw_adv = np.array([r["raw_adv"] for r in fold_rows], dtype=float)
    arr_raw_gap = np.array([r["raw_control_gap"] for r in fold_rows], dtype=float)

    arr_cal_auc = np.array([r["cal_auc"] for r in fold_rows], dtype=float)
    arr_cal_adv = np.array([r["cal_adv"] for r in fold_rows], dtype=float)
    arr_cal_gap = np.array([r["cal_control_gap"] for r in fold_rows], dtype=float)

    summary = {
        "n_rows": int(len(d)),
        "n_folds": n_folds,
        "raw": {
            "mean_auc": float(np.mean(arr_raw_auc)),
            "mean_adv": float(np.mean(arr_raw_adv)),
            "mean_control_gap": float(np.mean(arr_raw_gap)),
        },
        "calibrated": {
            "mean_auc": float(np.mean(arr_cal_auc)),
            "mean_adv": float(np.mean(arr_cal_adv)),
            "mean_control_gap": float(np.mean(arr_cal_gap)),
            "prob_adv_positive": float(np.mean(arr_cal_adv > 0.0)),
        },
        "improvement": {
            "delta_auc": float(np.mean(arr_cal_auc) - np.mean(arr_raw_auc)),
            "delta_adv": float(np.mean(arr_cal_adv) - np.mean(arr_raw_adv)),
            "gap_reduction": float(np.mean(arr_raw_gap) - np.mean(arr_cal_gap)),
            "gap_reduction_ratio": float((np.mean(arr_raw_gap) - np.mean(arr_cal_gap)) / max(np.mean(arr_raw_gap), 1e-12)),
        },
    }

    pass_auc = summary["calibrated"]["mean_auc"] >= 0.85
    pass_adv = summary["calibrated"]["mean_adv"] >= 0.50 and summary["calibrated"]["prob_adv_positive"] >= 0.90
    pass_gap = summary["calibrated"]["mean_control_gap"] <= 0.0010
    pass_gain = summary["improvement"]["gap_reduction_ratio"] >= 0.60

    if pass_auc and pass_adv and pass_gap and pass_gain:
        verdict = "GW_CONTROL_CALIBRATED_OBJECTIVE_SUPPORTED"
    elif pass_auc and pass_adv and pass_gain:
        verdict = "GW_CONTROL_CALIBRATED_OBJECTIVE_PARTIAL"
    else:
        verdict = "GW_CONTROL_CALIBRATED_OBJECTIVE_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "pass_flags": {
            "auc_support": bool(pass_auc),
            "adv_support": bool(pass_adv),
            "control_gap_support": bool(pass_gap),
            "gap_reduction_support": bool(pass_gain),
        },
        "verdict": verdict,
        "folds": fold_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1836: GW CONTROL-CALIBRATED OBJECTIVE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Rows/Folds: {summary['n_rows']}/{summary['n_folds']}",
        f"- Raw mean AUC/adv/gap: {summary['raw']['mean_auc']:.4f} / {summary['raw']['mean_adv']:.4f} / {summary['raw']['mean_control_gap']:.6f}",
        f"- Cal mean AUC/adv/gap: {summary['calibrated']['mean_auc']:.4f} / {summary['calibrated']['mean_adv']:.4f} / {summary['calibrated']['mean_control_gap']:.6f}",
        f"- Gap reduction ratio: {summary['improvement']['gap_reduction_ratio']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- auc_support: {pass_auc}",
        f"- adv_support: {pass_adv}",
        f"- control_gap_support: {pass_gap}",
        f"- gap_reduction_support: {pass_gain}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1836] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1836] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
