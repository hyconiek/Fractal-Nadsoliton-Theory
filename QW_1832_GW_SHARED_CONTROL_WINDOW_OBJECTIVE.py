#!/usr/bin/env python3
"""
QW-1832: GW shared-vs-control objective on event-windowed features.

Uses QW-1831 feature table to test whether H1-L1 windows are statistically
separable from control windows (H1-V1, L1-V1) with blocked CV.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
IN_CSV = ROOT / "gw1831_window_features.csv"
OUT_JSON = ROOT / "report_qw1832_gw_shared_control_window_objective.json"
OUT_MD = ROOT / "RAPORT_QW1832_GW_SHARED_CONTROL_WINDOW_OBJECTIVE.md"


def auc_from_scores(y: np.ndarray, s: np.ndarray) -> float:
    """AUC via rank formulation (Mann-Whitney)."""
    y = y.astype(int)
    n1 = int(np.sum(y == 1))
    n0 = int(np.sum(y == 0))
    if n1 == 0 or n0 == 0:
        return np.nan
    order = np.argsort(s)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    rank_sum_pos = float(np.sum(ranks[y == 1]))
    auc = (rank_sum_pos - n1 * (n1 + 1) / 2.0) / (n1 * n0)
    return float(auc)


def standardize(train_x: np.ndarray, test_x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mu = np.mean(train_x, axis=0, keepdims=True)
    sd = np.std(train_x, axis=0, keepdims=True)
    sd = np.where(sd <= 1e-12, 1.0, sd)
    return (train_x - mu) / sd, (test_x - mu) / sd


def fisher_linear_discriminant(x: np.ndarray, y: np.ndarray, reg: float = 1e-2) -> Tuple[np.ndarray, float]:
    x1 = x[y == 1]
    x0 = x[y == 0]
    m1 = np.mean(x1, axis=0)
    m0 = np.mean(x0, axis=0)
    s1 = np.cov(x1, rowvar=False)
    s0 = np.cov(x0, rowvar=False)
    if np.ndim(s1) == 0:
        s1 = np.array([[float(s1)]])
        s0 = np.array([[float(s0)]])
    sw = s1 + s0 + reg * np.eye(s1.shape[0])
    w = np.linalg.solve(sw, (m1 - m0))
    b = -0.5 * float(np.dot(w, (m1 + m0)))
    return w, b


def balanced_accuracy(y: np.ndarray, yhat: np.ndarray) -> float:
    y = y.astype(int)
    yhat = yhat.astype(int)
    tp = np.sum((y == 1) & (yhat == 1))
    fn = np.sum((y == 1) & (yhat == 0))
    tn = np.sum((y == 0) & (yhat == 0))
    fp = np.sum((y == 0) & (yhat == 1))
    tpr = tp / max(tp + fn, 1)
    tnr = tn / max(tn + fp, 1)
    return float(0.5 * (tpr + tnr))


def main() -> None:
    if not IN_CSV.exists():
        raise FileNotFoundError(f"Missing input CSV: {IN_CSV}")

    df = pd.read_csv(IN_CSV)
    req_cols = ["pair", "window_idx", "max_abs_corr", "corr_at_0ms", "corr_at_10ms", "mean_abs_corr", "best_lag_ms"]
    for c in req_cols:
        if c not in df.columns:
            raise RuntimeError(f"Missing column: {c}")

    data = df.copy()
    data["abs_best_lag_ms"] = np.abs(data["best_lag_ms"].astype(float))
    data["label_shared"] = (data["pair"] == "H1-L1").astype(int)

    feat_cols = ["max_abs_corr", "corr_at_0ms", "corr_at_10ms", "mean_abs_corr", "abs_best_lag_ms"]

    # Drop any residual NaN rows (should be none after QW-1831 fix)
    data = data.dropna(subset=feat_cols + ["window_idx", "label_shared"])

    n_folds = 5
    fold_rows: List[Dict[str, float]] = []

    for fold in range(n_folds):
        test_mask = (data["window_idx"].astype(int) % n_folds) == fold
        train_df = data[~test_mask]
        test_df = data[test_mask]

        x_train = train_df[feat_cols].to_numpy(dtype=float)
        y_train = train_df["label_shared"].to_numpy(dtype=int)
        x_test = test_df[feat_cols].to_numpy(dtype=float)
        y_test = test_df["label_shared"].to_numpy(dtype=int)

        x_train_s, x_test_s = standardize(x_train, x_test)
        w, b = fisher_linear_discriminant(x_train_s, y_train, reg=5e-3)

        s_train = x_train_s @ w + b
        s_test = x_test_s @ w + b

        auc = auc_from_scores(y_test, s_test)
        yhat = (s_test > 0.0).astype(int)
        bacc = balanced_accuracy(y_test, yhat)

        # Shared-vs-control prevalence advantage at high-score tail
        thr = float(np.quantile(s_train[y_train == 0], 0.90))
        p_shared = float(np.mean(s_test[y_test == 1] > thr))
        p_control = float(np.mean(s_test[y_test == 0] > thr))
        adv = p_shared - p_control

        fold_rows.append(
            {
                "fold": fold,
                "n_train": int(len(train_df)),
                "n_test": int(len(test_df)),
                "auc": float(auc),
                "balanced_accuracy": float(bacc),
                "thr_control_q90": thr,
                "p_shared_above_thr": p_shared,
                "p_control_above_thr": p_control,
                "prevalence_advantage": float(adv),
            }
        )

    arr_auc = np.array([r["auc"] for r in fold_rows], dtype=float)
    arr_ba = np.array([r["balanced_accuracy"] for r in fold_rows], dtype=float)
    arr_adv = np.array([r["prevalence_advantage"] for r in fold_rows], dtype=float)

    summary = {
        "n_rows": int(len(data)),
        "n_folds": n_folds,
        "feature_columns": feat_cols,
        "mean_auc": float(np.mean(arr_auc)),
        "std_auc": float(np.std(arr_auc)),
        "mean_balanced_accuracy": float(np.mean(arr_ba)),
        "std_balanced_accuracy": float(np.std(arr_ba)),
        "mean_prevalence_advantage": float(np.mean(arr_adv)),
        "std_prevalence_advantage": float(np.std(arr_adv)),
        "prob_adv_positive": float(np.mean(arr_adv > 0.0)),
    }

    pass_auc = summary["mean_auc"] >= 0.70
    pass_ba = summary["mean_balanced_accuracy"] >= 0.65
    pass_adv = summary["mean_prevalence_advantage"] >= 0.20 and summary["prob_adv_positive"] >= 0.80

    if pass_auc and pass_ba and pass_adv:
        verdict = "GW_WINDOW_OBJECTIVE_SUPPORTED"
    elif pass_auc and pass_ba:
        verdict = "GW_WINDOW_OBJECTIVE_PARTIAL"
    else:
        verdict = "GW_WINDOW_OBJECTIVE_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "pass_flags": {
            "auc_support": bool(pass_auc),
            "balanced_accuracy_support": bool(pass_ba),
            "prevalence_advantage_support": bool(pass_adv),
        },
        "verdict": verdict,
        "folds": fold_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1832: GW SHARED-CONTROL WINDOW OBJECTIVE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Rows: {summary['n_rows']}",
        f"- Folds: {summary['n_folds']}",
        f"- Mean AUC: {summary['mean_auc']:.4f} +/- {summary['std_auc']:.4f}",
        f"- Mean balanced accuracy: {summary['mean_balanced_accuracy']:.4f} +/- {summary['std_balanced_accuracy']:.4f}",
        f"- Mean prevalence advantage: {summary['mean_prevalence_advantage']:.4f} +/- {summary['std_prevalence_advantage']:.4f}",
        f"- P(adv>0): {summary['prob_adv_positive']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- auc_support: {pass_auc}",
        f"- balanced_accuracy_support: {pass_ba}",
        f"- prevalence_advantage_support: {pass_adv}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1832] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1832] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
