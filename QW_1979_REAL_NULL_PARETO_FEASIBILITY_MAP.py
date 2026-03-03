#!/usr/bin/env python3
"""
QW-1979: Real-vs-null Pareto feasibility map.

Purpose:
- quantify whether current operator class can satisfy closure-like constraints
  on both real pass-rate and null robustness simultaneously.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1978 = ROOT / "report_qw1978_worstcase_null_contrast_frontier.json"
OUT_JSON = ROOT / "report_qw1979_real_null_pareto_feasibility_map.json"
OUT_MD = ROOT / "RAPORT_QW1979_REAL_NULL_PARETO_FEASIBILITY_MAP.md"


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


def gw_metrics(s_hl: np.ndarray, s_hv: np.ndarray, s_lv: np.ndarray) -> Dict[str, float]:
    s_ctrl = np.concatenate([s_hv, s_lv])
    q90 = float(np.quantile(s_ctrl, 0.90))
    auc = rank_auc_pos_gt_neg(s_hl, s_ctrl)
    adv = float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90))
    sep = float(np.median(s_hl) - np.median(s_ctrl))
    gap = float(abs(np.median(s_hv) - np.median(s_lv)))
    return {
        "auc_h1l1_vs_ctrl": auc,
        "adv_shared_minus_ctrl_q90": adv,
        "sep_median_h1l1_minus_ctrl": sep,
        "control_median_gap": gap,
    }


def gw_flags(g: Dict[str, float], thr: Dict[str, float]) -> Dict[str, bool]:
    return {
        "gw_sep_ge_min": bool(g["sep_median_h1l1_minus_ctrl"] >= thr["gw_sep_min"]),
        "gw_adv_ge_min": bool(g["adv_shared_minus_ctrl_q90"] >= thr["gw_adv_min"]),
        "gw_auc_ge_min": bool(g["auc_h1l1_vs_ctrl"] >= thr["gw_auc_min"]),
        "gw_control_gap_le_max": bool(g["control_median_gap"] <= thr["gw_control_gap_max"]),
    }


def bootstrap_pass_rate(
    s_hl: np.ndarray,
    s_hv: np.ndarray,
    s_lv: np.ndarray,
    thr: Dict[str, float],
    n_boot: int,
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    n_hl, n_hv, n_lv = len(s_hl), len(s_hv), len(s_lv)
    pass_count = 0
    for _ in range(n_boot):
        b_hl = s_hl[rng.integers(0, n_hl, size=n_hl, endpoint=False)]
        b_hv = s_hv[rng.integers(0, n_hv, size=n_hv, endpoint=False)]
        b_lv = s_lv[rng.integers(0, n_lv, size=n_lv, endpoint=False)]
        if all(gw_flags(gw_metrics(b_hl, b_hv, b_lv), thr).values()):
            pass_count += 1
    return float(pass_count / n_boot)


def pareto_front(rows: List[Dict[str, float]]) -> List[Dict[str, float]]:
    # Maximize real_rate, minimize null_p90.
    pts = sorted(rows, key=lambda r: (-r["real_rate"], r["null_p90"]))
    front: List[Dict[str, float]] = []
    best_null = float("inf")
    for r in pts:
        if r["null_p90"] < best_null:
            front.append(r)
            best_null = r["null_p90"]
    return front


def verdict(strict_count: int, medium_count: int, relaxed_count: int) -> Tuple[str, str]:
    if strict_count > 0:
        return (
            "PARETO_STRICT_FEASIBLE",
            "PROMOTE_STRICT_CANDIDATES_TO_FULL_GATE",
        )
    if medium_count > 0:
        return (
            "PARETO_MEDIUM_FEASIBLE_ONLY",
            "NEED_BETTER_NULL_SUPPRESSION_FOR_STRICT_CLOSURE",
        )
    if relaxed_count > 0:
        return (
            "PARETO_RELAXED_ONLY",
            "CURRENT_OPERATOR_CLASS_INSUFFICIENT_FOR_CLOSURE",
        )
    return (
        "PARETO_NO_FEASIBLE_REGION",
        "REWORK_OPERATOR_CLASS_REQUIRED",
    )


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1978 = json.loads(IN_QW1978.read_text(encoding="utf-8"))
    thr = json.loads((ROOT / "report_qw1969_bootstrap_robust_recenter_search.json").read_text(encoding="utf-8"))["thresholds"]

    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    center = r1978["best"]

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
    pairs = df["pair"].map(pair_map).to_numpy(dtype=int)
    feats = df[["max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]].to_numpy(dtype=float)

    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    kvec = np.cos(kernel["omega"] * d + kernel["phi"]) / (1.0 + kernel["beta"] * (d**kernel["eta"]))
    w_raw = (np.abs(kvec) ** params["p_amp"]) * (d**params["r_dist"])
    w = w_raw / np.sum(w_raw)
    base_score = feats @ w

    lag_s = df["best_lag_ms"].to_numpy(dtype=float) * 1e-3
    lag_phase_sin = np.sin(kernel["omega"] * lag_s + kernel["phi"])
    lag_phase_cos = np.cos(kernel["omega"] * lag_s + kernel["phi"])
    corr0 = df["corr_at_0ms"].to_numpy(dtype=float)
    corr10 = df["corr_at_10ms"].to_numpy(dtype=float)
    mean_abs = df["mean_abs_corr"].to_numpy(dtype=float)

    pair_sign = np.where(pairs == 1, 1.0, np.where(pairs == 2, -1.0, 0.0))
    ctrl_mask = np.where(pairs == 0, 0.0, 1.0)
    ctrl_idx = np.where(pairs != 0)[0]
    n_ctrl = len(ctrl_idx)
    n_plus = int(np.sum(pairs == 1))

    c1_raw = pair_sign * lag_phase_sin * (corr0 - corr10)
    c2_raw = ctrl_mask * lag_phase_cos * (corr0 + corr10 + mean_abs)
    c1 = c1_raw / max(float(np.std(c1_raw)), 1e-12)
    c2 = c2_raw / max(float(np.std(c2_raw)), 1e-12)

    # Random local map around QW-1978 best.
    rng = np.random.default_rng(90790)
    n_samples = 500
    xi1 = np.clip(rng.normal(float(center["xi1"]), 0.00018, size=n_samples), 0.0002, 0.0016)
    xi2 = np.clip(rng.normal(float(center["xi2"]), 0.00022, size=n_samples), -0.0013, -0.00001)

    rows = []
    rng_null = np.random.default_rng(91790)
    for i in range(n_samples):
        a = float(xi1[i])
        b = float(xi2[i])
        score = base_score + a * c1 + b * c2
        s_hl = score[pairs == 0]
        s_hv = score[pairs == 1]
        s_lv = score[pairs == 2]
        g = gw_metrics(s_hl, s_hv, s_lv)
        if not all(gw_flags(g, thr).values()):
            continue

        real_rate = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=80, seed=92000 + i)

        null_rates = []
        for j in range(4):
            signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
            rng_null.shuffle(signs)
            rand_sign = np.zeros_like(pair_sign)
            rand_sign[ctrl_idx] = signs
            c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
            c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
            score_n = base_score + a * c1n + b * c2
            rb = bootstrap_pass_rate(
                score_n[pairs == 0], score_n[pairs == 1], score_n[pairs == 2], thr, n_boot=60, seed=93000 + i * 10 + j
            )
            null_rates.append(float(rb))

        null_mean = float(np.mean(null_rates))
        null_p90 = float(np.quantile(np.array(null_rates, dtype=float), 0.90))
        rows.append(
            {
                "xi1": a,
                "xi2": b,
                "real_rate": real_rate,
                "null_mean": null_mean,
                "null_p90": null_p90,
                "contrast_mean": float(real_rate - null_mean),
                "contrast_conservative": float(real_rate - null_p90),
            }
        )

    front = pareto_front(rows)

    strict = [r for r in rows if r["real_rate"] >= 0.90 and r["null_p90"] <= 0.45]
    medium = [r for r in rows if r["real_rate"] >= 0.85 and r["null_p90"] <= 0.50]
    relaxed = [r for r in rows if r["real_rate"] >= 0.80 and r["null_p90"] <= 0.55]

    best_cons = max(rows, key=lambda r: r["contrast_conservative"]) if rows else None
    best_mean = max(rows, key=lambda r: r["contrast_mean"]) if rows else None

    v, nxt = verdict(len(strict), len(medium), len(relaxed))

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1978.name],
        "scan": {
            "n_samples": int(n_samples),
            "gw_threshold_pass_rows": int(len(rows)),
        },
        "feasibility_counts": {
            "strict_real_ge_0p90_nullp90_le_0p45": int(len(strict)),
            "medium_real_ge_0p85_nullp90_le_0p50": int(len(medium)),
            "relaxed_real_ge_0p80_nullp90_le_0p55": int(len(relaxed)),
        },
        "best_conservative_contrast": best_cons,
        "best_mean_contrast": best_mean,
        "pareto_front_top20": front[:20],
        "verdict": v,
        "required_next_step": nxt,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1979: REAL-NULL PARETO FEASIBILITY MAP",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Feasibility Counts",
        f"- strict (real>=0.90, null_p90<=0.45): {len(strict)}",
        f"- medium (real>=0.85, null_p90<=0.50): {len(medium)}",
        f"- relaxed (real>=0.80, null_p90<=0.55): {len(relaxed)}",
        "",
        "## Best Conservative Contrast",
        (
            f"- real/null_mean/null_p90/cons: "
            f"{100.0 * best_cons['real_rate']:.2f}%/"
            f"{100.0 * best_cons['null_mean']:.2f}%/"
            f"{100.0 * best_cons['null_p90']:.2f}%/"
            f"{100.0 * best_cons['contrast_conservative']:.2f} pp"
        )
        if best_cons
        else "- none",
        "",
        "## Required Next Step",
        f"- {nxt}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1979] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1979] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1979] verdict={v}")


if __name__ == "__main__":
    main()
