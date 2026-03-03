#!/usr/bin/env python3
"""
QW-1841: Permutation-null test for GW control-calibrated objective (QW-1836).

Goal:
- estimate how extreme the observed calibrated GW metrics are under a null
  where pair labels carry no information (label permutation test).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
IN_CSV = ROOT / "gw1831_window_features.csv"
OUT_JSON = ROOT / "report_qw1841_gw_calibrated_permutation_null.json"
OUT_MD = ROOT / "RAPORT_QW1841_GW_CALIBRATED_PERMUTATION_NULL.md"


def auc_pos_gt_neg(pos: np.ndarray, neg: np.ndarray) -> float:
    y = np.concatenate([np.ones(len(pos), dtype=int), np.zeros(len(neg), dtype=int)])
    s = np.concatenate([pos, neg])
    n1 = len(pos)
    n0 = len(neg)
    if n1 == 0 or n0 == 0:
        return float("nan")
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
    if len(sh) == 0 or len(hv) == 0 or len(lv) == 0:
        return {
            "auc": float("nan"),
            "adv": float("nan"),
            "gap": float("nan"),
        }

    ctrl = np.concatenate([hv, lv])
    auc = auc_pos_gt_neg(sh, ctrl)
    q90 = float(np.quantile(ctrl, 0.90))
    p_shared = float(np.mean(sh > q90))
    p_ctrl = float(np.mean(ctrl > q90))
    adv = float(p_shared - p_ctrl)
    gap = float(abs(np.median(hv) - np.median(lv)))
    return {"auc": float(auc), "adv": adv, "gap": gap}


def evaluate_calibrated(pair_labels: np.ndarray, window_idx: np.ndarray, score: np.ndarray, n_folds: int = 5) -> Tuple[float, float, float]:
    aucs = []
    advs = []
    gaps = []

    for fold in range(n_folds):
        tr = (window_idx % n_folds) != fold
        te = ~tr

        p_tr = pair_labels[tr]
        s_tr = score[tr]
        p_te = pair_labels[te]
        s_te = score[te].copy()

        hv_tr = s_tr[p_tr == "H1-V1"]
        lv_tr = s_tr[p_tr == "L1-V1"]
        c_tr = s_tr[p_tr != "H1-L1"]
        if len(hv_tr) == 0 or len(lv_tr) == 0 or len(c_tr) == 0:
            return (float("nan"), float("nan"), float("nan"))

        off_hv = float(np.median(hv_tr) - np.median(c_tr))
        off_lv = float(np.median(lv_tr) - np.median(c_tr))

        s_te[p_te == "H1-V1"] = s_te[p_te == "H1-V1"] - off_hv
        s_te[p_te == "L1-V1"] = s_te[p_te == "L1-V1"] - off_lv

        m = metrics_from_scores(p_te, s_te)
        aucs.append(m["auc"])
        advs.append(m["adv"])
        gaps.append(m["gap"])

    return float(np.mean(aucs)), float(np.mean(advs)), float(np.mean(gaps))


def empirical_p_ge(null: np.ndarray, obs: float) -> float:
    return float((1 + np.sum(null >= obs)) / (len(null) + 1))


def empirical_p_le(null: np.ndarray, obs: float) -> float:
    return float((1 + np.sum(null <= obs)) / (len(null) + 1))


def main() -> None:
    if not IN_CSV.exists():
        raise FileNotFoundError(f"Missing input CSV: {IN_CSV}")

    df = pd.read_csv(IN_CSV)
    req = ["pair", "window_idx", "max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]
    for c in req:
        if c not in df.columns:
            raise RuntimeError(f"Missing column: {c}")

    d = df.dropna(subset=req).copy()
    pair = d["pair"].to_numpy()
    window_idx = d["window_idx"].to_numpy(dtype=int)
    s_raw = score_raw(d)

    obs_auc, obs_adv, obs_gap = evaluate_calibrated(pair, window_idx, s_raw, n_folds=5)

    n_perm = 5000
    rng = np.random.default_rng(21926)

    null_auc = np.empty(n_perm, dtype=float)
    null_adv = np.empty(n_perm, dtype=float)
    null_gap = np.empty(n_perm, dtype=float)

    for i in range(n_perm):
        p_perm = rng.permutation(pair)
        a, v, g = evaluate_calibrated(p_perm, window_idx, s_raw, n_folds=5)
        null_auc[i] = a
        null_adv[i] = v
        null_gap[i] = g

    p_auc = empirical_p_ge(null_auc, obs_auc)
    p_adv = empirical_p_ge(null_adv, obs_adv)
    p_gap = empirical_p_le(null_gap, obs_gap)

    pass_auc = p_auc <= 0.01
    pass_adv = p_adv <= 0.01
    pass_gap = p_gap <= 0.01

    if pass_auc and pass_adv and pass_gap:
        verdict = "GW_CALIBRATED_OBJECTIVE_NULL_REJECTED"
    elif pass_auc and pass_adv:
        verdict = "GW_CALIBRATED_OBJECTIVE_PARTIAL_NULL_REJECTION"
    else:
        verdict = "GW_CALIBRATED_OBJECTIVE_NULL_NOT_REJECTED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_rows": int(len(d)),
        "n_perm": n_perm,
        "observed": {
            "mean_auc": obs_auc,
            "mean_adv": obs_adv,
            "mean_control_gap": obs_gap,
        },
        "null": {
            "auc": {
                "mean": float(np.mean(null_auc)),
                "std": float(np.std(null_auc)),
                "q95": float(np.quantile(null_auc, 0.95)),
            },
            "adv": {
                "mean": float(np.mean(null_adv)),
                "std": float(np.std(null_adv)),
                "q95": float(np.quantile(null_adv, 0.95)),
            },
            "gap": {
                "mean": float(np.mean(null_gap)),
                "std": float(np.std(null_gap)),
                "q05": float(np.quantile(null_gap, 0.05)),
            },
        },
        "p_values": {
            "auc_right_tail": p_auc,
            "adv_right_tail": p_adv,
            "gap_left_tail": p_gap,
        },
        "pass_flags": {
            "auc_null_rejected": bool(pass_auc),
            "adv_null_rejected": bool(pass_adv),
            "gap_null_rejected": bool(pass_gap),
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1841: GW CALIBRATED PERMUTATION NULL",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Rows/Permutations: {out['n_rows']}/{n_perm}",
        f"- Obserwowane AUC/ADV/GAP: {obs_auc:.4f} / {obs_adv:.4f} / {obs_gap:.6f}",
        (
            "- Null means AUC/ADV/GAP: "
            f"{out['null']['auc']['mean']:.4f} / {out['null']['adv']['mean']:.4f} / {out['null']['gap']['mean']:.6f}"
        ),
        (
            "- p-values (AUC right, ADV right, GAP left): "
            f"{p_auc:.6f}, {p_adv:.6f}, {p_gap:.6f}"
        ),
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags (alpha=0.01)",
        f"- auc_null_rejected: {pass_auc}",
        f"- adv_null_rejected: {pass_adv}",
        f"- gap_null_rejected: {pass_gap}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1841] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1841] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
