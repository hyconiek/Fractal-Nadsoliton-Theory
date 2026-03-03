#!/usr/bin/env python3
"""
QW-1983: Fold-robust operator search (single operator, no fold retune).
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
IN_QW1982 = ROOT / "report_qw1982_temporal_fold_transfer_audit.json"
OUT_JSON = ROOT / "report_qw1983_fold_robust_operator_search.json"
OUT_MD = ROOT / "RAPORT_QW1983_FOLD_ROBUST_OPERATOR_SEARCH.md"

# Compute budget (balanced for turnaround + statistical stability).
N_SEARCH = 140
FAST_REAL_BOOT = 120
FAST_NULL_TRIALS = 8
FAST_NULL_BOOT = 40
FULL_REAL_BOOT = 600
SHORTLIST_SIZE = 8


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


def build_fold_channels(
    df_fold: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
    xi1: float,
    xi3: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
    pairs = df_fold["pair"].map(pair_map).to_numpy(dtype=int)

    feats = df_fold[["max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]].to_numpy(dtype=float)
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    kvec = np.cos(kernel["omega"] * d + kernel["phi"]) / (1.0 + kernel["beta"] * (d**kernel["eta"]))
    w_raw = (np.abs(kvec) ** params["p_amp"]) * (d**params["r_dist"])
    w = w_raw / np.sum(w_raw)
    base_score = feats @ w

    lag_s = df_fold["best_lag_ms"].to_numpy(dtype=float) * 1e-3
    lag_phase_sin = np.sin(kernel["omega"] * lag_s + kernel["phi"])
    lag_phase_cos = np.cos(kernel["omega"] * lag_s + kernel["phi"])
    corr0 = df_fold["corr_at_0ms"].to_numpy(dtype=float)
    corr10 = df_fold["corr_at_10ms"].to_numpy(dtype=float)
    mean_abs = df_fold["mean_abs_corr"].to_numpy(dtype=float)

    pair_sign = np.where(pairs == 1, 1.0, np.where(pairs == 2, -1.0, 0.0))
    c1_raw = pair_sign * lag_phase_sin * (corr0 - corr10)
    c3_raw = pair_sign * lag_phase_cos * (mean_abs + np.abs(corr0) + np.abs(corr10))
    c1 = c1_raw / max(float(np.std(c1_raw)), 1e-12)
    c3 = c3_raw / max(float(np.std(c3_raw)), 1e-12)

    s = base_score + xi1 * c1 + xi3 * c3
    return s[pairs == 0], s[pairs == 1], s[pairs == 2], pairs, c1, c3


def fold_null_stats_fast(
    s_base: np.ndarray,
    pairs: np.ndarray,
    c1: np.ndarray,
    c3: np.ndarray,
    xi1: float,
    xi3: float,
    thr: Dict[str, float],
    seed: int,
) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    ctrl_idx = np.where(pairs != 0)[0]
    n_ctrl = len(ctrl_idx)
    n_plus = int(np.sum(pairs == 1))
    pair_sign = np.where(pairs == 1, 1.0, np.where(pairs == 2, -1.0, 0.0))

    random_rates = []
    for i in range(FAST_NULL_TRIALS):
        signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
        rng.shuffle(signs)
        rand_sign = np.zeros_like(pair_sign)
        rand_sign[ctrl_idx] = signs

        c1n_raw = rand_sign * c1
        c3n_raw = rand_sign * c3
        c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
        c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
        s = s_base + xi1 * c1n + xi3 * c3n
        s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]
        rb = bootstrap_pass_rate(
            s_hl,
            s_hv,
            s_lv,
            thr,
            n_boot=FAST_NULL_BOOT,
            seed=seed + 100 + i,
        )
        random_rates.append(float(rb))

    return float(np.mean(random_rates)), float(np.quantile(np.array(random_rates, dtype=float), 0.90))


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1982 = json.loads(IN_QW1982.read_text(encoding="utf-8"))
    thr = json.loads((ROOT / "report_qw1969_bootstrap_robust_recenter_search.json").read_text(encoding="utf-8"))["thresholds"]
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]

    # center around failed candidate from QW-1982
    xi1_center = float(r1982["candidate"]["xi1"])
    xi3_center = float(r1982["candidate"]["xi3"])

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    rng = np.random.default_rng(130000)
    n_search = N_SEARCH
    xi1_samples = np.clip(rng.normal(xi1_center, 0.00035, size=n_search), 0.0001, 0.0025)
    xi3_samples = np.clip(rng.normal(xi3_center, 0.00035, size=n_search), -0.0002, 0.0022)

    candidates = []
    for i in range(n_search):
        xi1 = float(xi1_samples[i])
        xi3 = float(xi3_samples[i])
        fold_rows = []
        for f, dff in enumerate(fold_dfs):
            s_hl, s_hv, s_lv, pairs, c1, c3 = build_fold_channels(dff, kernel, params, xi1, xi3)
            g = gw_metrics(s_hl, s_hv, s_lv)
            flags = gw_flags(g, thr)
            real_fast = bootstrap_pass_rate(
                s_hl,
                s_hv,
                s_lv,
                thr,
                n_boot=FAST_REAL_BOOT,
                seed=131000 + i * 10 + f,
            )
            # reconstruct base signal for null generation
            s_base = s_hl.copy()  # placeholder not used directly
            # we need full s_base vector; derive from parts:
            # simpler: rebuild from channels by insertion is cumbersome, so recompute:
            pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
            pairs_vec = dff["pair"].map(pair_map).to_numpy(dtype=int)
            # c1 and c3 are full vectors already from build_fold_channels last return
            # reconstruct s from channels:
            s_full = np.zeros(len(pairs_vec), dtype=float)
            s_full[pairs_vec == 0] = s_hl
            s_full[pairs_vec == 1] = s_hv
            s_full[pairs_vec == 2] = s_lv
            null_mean, null_p90 = fold_null_stats_fast(
                s_base=s_full - xi1 * c1 - xi3 * c3,
                pairs=pairs_vec,
                c1=c1,
                c3=c3,
                xi1=xi1,
                xi3=xi3,
                thr=thr,
                seed=132000 + i * 10 + f,
            )
            fold_rows.append(
                {
                    "fold": f,
                    "real_fast": real_fast,
                    "null_mean_fast": null_mean,
                    "null_p90_fast": null_p90,
                    "det_flags": flags,
                }
            )

        min_real = float(min(r["real_fast"] for r in fold_rows))
        max_null_p90 = float(max(r["null_p90_fast"] for r in fold_rows))
        robust_score = float(min_real - 0.9 * max_null_p90)
        candidates.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "fold_fast": fold_rows,
                "min_real_fast": min_real,
                "max_null_p90_fast": max_null_p90,
                "robust_score": robust_score,
            }
        )
        if (i + 1) % 20 == 0:
            print(f"[QW-1983] search progress: {i + 1}/{n_search}", flush=True)

    candidates.sort(key=lambda x: x["robust_score"], reverse=True)
    shortlist = candidates[:SHORTLIST_SIZE]
    print(f"[QW-1983] shortlist size: {len(shortlist)}", flush=True)

    # Full audit for shortlist
    final_rows = []
    for i, cand in enumerate(shortlist):
        xi1 = cand["xi1"]
        xi3 = cand["xi3"]
        full_folds = []
        for f, dff in enumerate(fold_dfs):
            s_hl, s_hv, s_lv, pairs, c1, c3 = build_fold_channels(dff, kernel, params, xi1, xi3)
            g = gw_metrics(s_hl, s_hv, s_lv)
            flags = gw_flags(g, thr)
            real_full = bootstrap_pass_rate(
                s_hl,
                s_hv,
                s_lv,
                thr,
                n_boot=FULL_REAL_BOOT,
                seed=133000 + i * 10 + f,
            )

            pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
            pairs_vec = dff["pair"].map(pair_map).to_numpy(dtype=int)
            s_full = np.zeros(len(pairs_vec), dtype=float)
            s_full[pairs_vec == 0] = s_hl
            s_full[pairs_vec == 1] = s_hv
            s_full[pairs_vec == 2] = s_lv
            null_mean, null_p90 = fold_null_stats_fast(
                s_base=s_full - xi1 * c1 - xi3 * c3,
                pairs=pairs_vec,
                c1=c1,
                c3=c3,
                xi1=xi1,
                xi3=xi3,
                thr=thr,
                seed=134000 + i * 10 + f,
            )

            fold_pass = bool(real_full >= 0.85 and null_p90 <= 0.40)
            full_folds.append(
                {
                    "fold": f,
                    "real_full": real_full,
                    "null_mean_full": null_mean,
                    "null_p90_full": null_p90,
                    "det_flags": flags,
                    "fold_pass": fold_pass,
                }
            )

        pass_count = int(sum(int(r["fold_pass"]) for r in full_folds))
        min_real = float(min(r["real_full"] for r in full_folds))
        max_null_p90 = float(max(r["null_p90_full"] for r in full_folds))
        final_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "pass_count": pass_count,
                "min_real_full": min_real,
                "max_null_p90_full": max_null_p90,
                "robust_score_full": float(min_real - 0.9 * max_null_p90),
                "fold_results": full_folds,
            }
        )

    final_rows.sort(key=lambda x: (x["pass_count"], x["robust_score_full"]), reverse=True)
    best = final_rows[0]

    verdict = (
        "FOLD_ROBUST_SEARCH_PASS_STRONG"
        if best["pass_count"] >= 4
        else "FOLD_ROBUST_SEARCH_PASS_PARTIAL"
        if best["pass_count"] >= 3
        else "FOLD_ROBUST_SEARCH_FAIL"
    )
    required = (
        "PROMOTE_BEST_TO_EXTERNAL_CONFIRMATORY_PREP"
        if best["pass_count"] >= 4
        else "INCREASE_FOLD_ROBUSTNESS_AND_REPEAT"
        if best["pass_count"] >= 3
        else "REWORK_OPERATOR_CLASS_OR_OBJECTIVE"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1982.name],
        "search": {"n_search": n_search, "shortlist_size": len(shortlist), "final_count": len(final_rows)},
        "best": best,
        "top8": final_rows,
        "verdict": verdict,
        "required_next_step": required,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1983: FOLD ROBUST OPERATOR SEARCH",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Best Candidate",
        f"- xi1/xi3: {best['xi1']:.6f}/{best['xi3']:.6f}",
        f"- pass_count: {best['pass_count']}/5",
        f"- min_real_full: {100.0 * best['min_real_full']:.2f}%",
        f"- max_null_p90_full: {100.0 * best['max_null_p90_full']:.2f}%",
        "",
        "## Required Next Step",
        f"- {required}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1983] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1983] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1983] verdict={verdict} pass_count={best['pass_count']}/5")


if __name__ == "__main__":
    main()
