#!/usr/bin/env python3
"""
QW-1982: Temporal fold transfer audit for QW-1981 candidate.

Goal:
- test transfer without retuning across temporal folds,
- quantify real-vs-null stability in each fold separately.
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
IN_QW1981 = ROOT / "report_qw1981_signed_dual_stress_and_leakage_audit.json"
OUT_JSON = ROOT / "report_qw1982_temporal_fold_transfer_audit.json"
OUT_MD = ROOT / "RAPORT_QW1982_TEMPORAL_FOLD_TRANSFER_AUDIT.md"


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


def build_scores(
    df_fold: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
    xi1: float,
    xi3: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    return s[pairs == 0], s[pairs == 1], s[pairs == 2]


def fold_null_stats(
    df_fold: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
    xi1: float,
    xi3: float,
    thr: Dict[str, float],
    seed: int,
) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
    pairs = df_fold["pair"].map(pair_map).to_numpy(dtype=int)
    ctrl_idx = np.where(pairs != 0)[0]
    n_ctrl = len(ctrl_idx)
    n_plus = int(np.sum(pairs == 1))

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

    random_rates = []
    det_pass = 0
    for i in range(40):
        signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
        rng.shuffle(signs)
        rand_sign = np.zeros_like(pair_sign)
        rand_sign[ctrl_idx] = signs

        c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
        c3n_raw = rand_sign * lag_phase_cos * (mean_abs + np.abs(corr0) + np.abs(corr10))
        c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
        c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
        s = base_score + xi1 * c1n + xi3 * c3n
        s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]
        g = gw_metrics(s_hl, s_hv, s_lv)
        if all(gw_flags(g, thr).values()):
            det_pass += 1
        rb = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=200, seed=seed + 1000 + i)
        random_rates.append(float(rb))

    adv_rates = []
    for i in range(12):
        best_candidate = None
        for _ in range(20):
            signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
            rng.shuffle(signs)
            rand_sign = np.zeros_like(pair_sign)
            rand_sign[ctrl_idx] = signs
            c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
            c3n_raw = rand_sign * lag_phase_cos * (mean_abs + np.abs(corr0) + np.abs(corr10))
            c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
            c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
            s = base_score + xi1 * c1n + xi3 * c3n
            s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]
            g = gw_metrics(s_hl, s_hv, s_lv)
            det_score = (
                1.2 * (g["auc_h1l1_vs_ctrl"] - thr["gw_auc_min"])
                + 1.0 * (g["adv_shared_minus_ctrl_q90"] - thr["gw_adv_min"])
                + 120.0 * (g["sep_median_h1l1_minus_ctrl"] - thr["gw_sep_min"])
                + 180.0 * (thr["gw_control_gap_max"] - g["control_median_gap"])
            )
            if (best_candidate is None) or (det_score > best_candidate["det_score"]):
                best_candidate = {"s_hl": s_hl, "s_hv": s_hv, "s_lv": s_lv, "det_score": det_score}
        rb = bootstrap_pass_rate(best_candidate["s_hl"], best_candidate["s_hv"], best_candidate["s_lv"], thr, n_boot=200, seed=seed + 2000 + i)
        adv_rates.append(float(rb))

    return {
        "null_det_pass_rate": float(det_pass / 40.0),
        "null_random_mean": float(np.mean(random_rates)),
        "null_random_p90": float(np.quantile(np.array(random_rates, dtype=float), 0.90)),
        "null_adversarial_mean": float(np.mean(adv_rates)),
        "null_adversarial_p90": float(np.quantile(np.array(adv_rates, dtype=float), 0.90)),
    }


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1981 = json.loads(IN_QW1981.read_text(encoding="utf-8"))
    thr = json.loads((ROOT / "report_qw1969_bootstrap_robust_recenter_search.json").read_text(encoding="utf-8"))["thresholds"]

    xi1 = float(r1981["candidate"]["xi1"])
    xi3 = float(r1981["candidate"]["xi3"])
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    # temporal stratification by window index (same fold per pair index position)
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5

    fold_rows = []
    for fold in range(5):
        df_fold = df[df["fold"] == fold].reset_index(drop=True)
        s_hl, s_hv, s_lv = build_scores(df_fold, kernel, params, xi1, xi3)
        g = gw_metrics(s_hl, s_hv, s_lv)
        flags = gw_flags(g, thr)
        real_boot = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=1200, seed=120000 + fold)
        null_stats = fold_null_stats(df_fold, kernel, params, xi1, xi3, thr, seed=121000 + fold * 100)

        fold_pass = bool(
            (real_boot >= 0.85)
            and (null_stats["null_random_p90"] <= 0.40)
            and (null_stats["null_adversarial_mean"] <= 0.45)
        )

        fold_rows.append(
            {
                "fold": int(fold),
                "n_rows": int(len(df_fold)),
                "deterministic_gw": g,
                "deterministic_flags": flags,
                "real_boot_pass_rate_1200": real_boot,
                **null_stats,
                "fold_pass": fold_pass,
            }
        )

    pass_count = sum(int(r["fold_pass"]) for r in fold_rows)
    verdict = (
        "TEMPORAL_TRANSFER_PASS_STRONG"
        if pass_count >= 4
        else "TEMPORAL_TRANSFER_PASS_PARTIAL"
        if pass_count >= 3
        else "TEMPORAL_TRANSFER_FAIL"
    )
    required = (
        "PREPARE_TRUE_EXTERNAL_CONFIRMATORY_EXECUTION"
        if pass_count >= 4
        else "RUN_ADDITIONAL_TRANSFER_HOLDOUTS_BEFORE_EXTERNAL"
        if pass_count >= 3
        else "REJECT_FREEZE_AND_REWORK_OPERATOR"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1981.name],
        "candidate": {"xi1": xi1, "xi3": xi3},
        "thresholds_transfer_gate": {
            "real_boot_min": 0.85,
            "null_random_p90_max": 0.40,
            "null_adversarial_mean_max": 0.45,
            "required_folds_to_pass": 4,
        },
        "fold_results": fold_rows,
        "summary": {
            "n_folds": 5,
            "pass_count": int(pass_count),
            "real_boot_mean": float(np.mean([r["real_boot_pass_rate_1200"] for r in fold_rows])),
            "null_random_p90_mean": float(np.mean([r["null_random_p90"] for r in fold_rows])),
            "null_adversarial_mean_mean": float(np.mean([r["null_adversarial_mean"] for r in fold_rows])),
        },
        "verdict": verdict,
        "required_next_step": required,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1982: TEMPORAL FOLD TRANSFER AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Candidate",
        f"- xi1/xi3: {xi1:.6f}/{xi3:.6f}",
        "",
        "## Fold Summary",
        f"- pass_count: {pass_count}/5",
        f"- mean real_boot: {100.0 * out['summary']['real_boot_mean']:.2f}%",
        f"- mean null_random_p90: {100.0 * out['summary']['null_random_p90_mean']:.2f}%",
        f"- mean null_adversarial_mean: {100.0 * out['summary']['null_adversarial_mean_mean']:.2f}%",
        "",
        "## Required Next Step",
        f"- {required}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1982] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1982] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1982] verdict={verdict} pass_count={pass_count}/5")


if __name__ == "__main__":
    main()
