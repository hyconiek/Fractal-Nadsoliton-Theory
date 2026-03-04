#!/usr/bin/env python3
"""
QW-2061: GW control-gap feasibility frontier (with locked flavor context).

Question:
- under current GW feature space and linear shared score, is strict control_gap <= 0.0025
  simultaneously feasible with auc/adv/sep thresholds?
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2061_gw_control_gap_feasibility_frontier.json"
OUT_MD = ROOT / "RAPORT_QW2061_GW_CONTROL_GAP_FEASIBILITY_FRONTIER.md"


def rank_auc_pos_gt_neg(pos: np.ndarray, neg: np.ndarray) -> float:
    y = np.concatenate([np.ones(len(pos), dtype=int), np.zeros(len(neg), dtype=int)])
    s = np.concatenate([pos, neg])
    n1 = len(pos)
    n0 = len(neg)
    order = np.argsort(s)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    rs = float(np.sum(ranks[y == 1]))
    return float((rs - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def gw_metrics(score: np.ndarray, pair: np.ndarray) -> Dict[str, float]:
    s_hl = score[pair == "H1-L1"]
    s_hv = score[pair == "H1-V1"]
    s_lv = score[pair == "L1-V1"]
    s_ctrl = np.concatenate([s_hv, s_lv])
    q90 = float(np.quantile(s_ctrl, 0.90))
    return {
        "auc_h1l1_vs_ctrl": float(rank_auc_pos_gt_neg(s_hl, s_ctrl)),
        "adv_shared_minus_ctrl_q90": float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90)),
        "sep_median_h1l1_minus_ctrl": float(np.median(s_hl) - np.median(s_ctrl)),
        "control_median_gap": float(abs(np.median(s_hv) - np.median(s_lv))),
    }


def simplex_weights(step: float = 0.05) -> List[Tuple[float, float, float, float]]:
    n = int(round(1.0 / step))
    out: List[Tuple[float, float, float, float]] = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            for k in range(n + 1 - i - j):
                l = n - i - j - k
                out.append((i * step, j * step, k * step, l * step))
    return out


def main() -> None:
    r2060 = json.loads((ROOT / "report_qw2060_locked_shared_flavor_basis_no_rescan_gate.json").read_text(encoding="utf-8"))
    flavor_locked = r2060["metrics"]["flavor"]

    thresholds = {
        "gw_sep_min": 0.0020,
        "gw_adv_min": 0.30,
        "gw_auc_min": 0.75,
        "gw_control_gap_max": 0.0025,
    }

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    f1 = df["max_abs_corr"].to_numpy(dtype=float)
    f2 = df["mean_abs_corr"].to_numpy(dtype=float)
    f3 = df["corr_at_0ms"].to_numpy(dtype=float)
    f4 = df["corr_at_10ms"].to_numpy(dtype=float)
    pair = df["pair"].astype(str).to_numpy()

    rows = []
    for w1, w2, w3, w4 in simplex_weights(step=0.05):
        score = w1 * f1 + w2 * f2 + w3 * f3 + w4 * f4
        m = gw_metrics(score, pair)
        flags = {
            "gw_sep_ge_min": bool(m["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
            "gw_adv_ge_min": bool(m["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
            "gw_auc_ge_min": bool(m["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
            "gw_control_gap_le_max": bool(m["control_median_gap"] <= thresholds["gw_control_gap_max"]),
        }
        row = {
            "weights": {
                "w_max_abs_corr": float(w1),
                "w_mean_abs_corr": float(w2),
                "w_corr_at_0ms": float(w3),
                "w_corr_at_10ms": float(w4),
            },
            "metrics": m,
            "flags": flags,
            "all_pass": bool(all(flags.values())),
            "soft_pass_no_gap": bool(flags["gw_sep_ge_min"] and flags["gw_adv_ge_min"] and flags["gw_auc_ge_min"]),
            "score_gap_primary": float(
                2000.0 * max(0.0, m["control_median_gap"] - thresholds["gw_control_gap_max"])
                + 100.0 * max(0.0, thresholds["gw_sep_min"] - m["sep_median_h1l1_minus_ctrl"])
                + 50.0 * max(0.0, thresholds["gw_adv_min"] - m["adv_shared_minus_ctrl_q90"])
                + 50.0 * max(0.0, thresholds["gw_auc_min"] - m["auc_h1l1_vs_ctrl"])
            ),
        }
        rows.append(row)

    rows_sorted = sorted(rows, key=lambda r: float(r["score_gap_primary"]))
    pass_count = int(sum(1 for r in rows if r["all_pass"]))
    soft_count = int(sum(1 for r in rows if r["soft_pass_no_gap"]))

    best = rows_sorted[0]
    best_gap = min(rows, key=lambda r: float(r["metrics"]["control_median_gap"]))
    best_soft = min(
        [r for r in rows if r["soft_pass_no_gap"]],
        key=lambda r: float(r["metrics"]["control_median_gap"]),
    )

    verdict = (
        "GW_CONTROL_GAP_FEASIBLE_IN_CURRENT_LINEAR_FEATURE_SPACE"
        if pass_count > 0
        else "GW_CONTROL_GAP_INFEASIBLE_UNDER_CURRENT_LINEAR_FEATURE_OPERATOR"
    )
    required_next = (
        "LOCK_BEST_FEASIBLE_GW_WEIGHTS_WITH_LOCKED_FLAVOR_CONTEXT"
        if pass_count > 0
        else "INTRODUCE_HIGHER_ORDER_OR_NONLINEAR_GW_OPERATOR_OR_REASSESS_GAP_THRESHOLD"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "locked_flavor_context_source": "report_qw2060_locked_shared_flavor_basis_no_rescan_gate.json",
        "locked_flavor_metrics": flavor_locked,
        "gw_thresholds": thresholds,
        "search": {
            "simplex_step": 0.05,
            "n_candidates": int(len(rows)),
            "pass_count_all": pass_count,
            "soft_count_no_gap": soft_count,
        },
        "best_row": best,
        "best_control_gap_row_global": best_gap,
        "best_control_gap_row_with_soft_constraints": best_soft,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2061: GW CONTROL-GAP FEASIBILITY FRONTIER",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- n_candidates: {out['search']['n_candidates']}",
        f"- pass_count_all: {pass_count}",
        f"- soft_count_no_gap: {soft_count}",
        "",
        "## Locked Flavor Context (from QW-2060)",
        (
            f"- CKM/PMNS mean rel%: "
            f"{flavor_locked['ckm_mean_rel_pct']:.3f}/{flavor_locked['pmns_mean_rel_pct']:.3f}"
        ),
        "",
        "## Best Row (primary score)",
        f"- weights: {best['weights']}",
        (
            f"- GW auc/adv/sep/gap: "
            f"{best['metrics']['auc_h1l1_vs_ctrl']:.4f}/"
            f"{best['metrics']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{best['metrics']['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{best['metrics']['control_median_gap']:.6f}"
        ),
        f"- all_pass: {best['all_pass']}",
        f"- soft_pass_no_gap: {best['soft_pass_no_gap']}",
        "",
        "## Best Global Control-Gap Row",
        f"- weights: {best_gap['weights']}",
        (
            f"- GW auc/adv/sep/gap: "
            f"{best_gap['metrics']['auc_h1l1_vs_ctrl']:.4f}/"
            f"{best_gap['metrics']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{best_gap['metrics']['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{best_gap['metrics']['control_median_gap']:.6f}"
        ),
        "",
        "## Best Control-Gap Under Soft Constraints",
        f"- weights: {best_soft['weights']}",
        (
            f"- GW auc/adv/sep/gap: "
            f"{best_soft['metrics']['auc_h1l1_vs_ctrl']:.4f}/"
            f"{best_soft['metrics']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{best_soft['metrics']['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{best_soft['metrics']['control_median_gap']:.6f}"
        ),
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2061] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2061] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2061] verdict={verdict} pass_count={pass_count}")


if __name__ == "__main__":
    main()
