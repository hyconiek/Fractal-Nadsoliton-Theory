#!/usr/bin/env python3
"""
QW-1993: Global score-compression hard-lock search.
Adds one global monotonic compression parameter gamma_c in score readout.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import ROOT, bootstrap_pass_rate
from QW_1988_TRI_BASIS_HARD_LOCK_STRESS import bootstrap_pass_rate_block
from QW_1989_CONSTRAINED_ADVERSARIAL_AUDIT import optimal_sign_sequence_with_flip_budget


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1992 = ROOT / "report_qw1992_multi_fold_null_ceiling_balance_search.json"
OUT_JSON = ROOT / "report_qw1993_global_score_compression_hard_lock_search.json"
OUT_MD = ROOT / "RAPORT_QW1993_GLOBAL_SCORE_COMPRESSION_HARD_LOCK_SEARCH.md"

TARGET_REAL_IID_MIN = 0.95
TARGET_REAL_BLOCK_MIN = 0.90
TARGET_NULL_RANDOM_P90_MAX = 0.40

FAST_REAL_IID_BOOT = 120
FAST_NULL_TRIALS = 6
FAST_NULL_BOOT = 30

FULL_REAL_IID_BOOT = 1400
FULL_REAL_BLOCK_BOOT = 900
FULL_REAL_BLOCK_LEN = 10
FULL_NULL_TRIALS = 18
FULL_NULL_BOOT = 80
FULL_ADV_BOOT = 350
FULL_ADV_FLIP_BUDGET = 4

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


def build_components(
    dff: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
    pairs = dff["pair"].map(pair_map).to_numpy(dtype=int)

    feats = dff[["max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]].to_numpy(dtype=float)
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    kvec = np.cos(kernel["omega"] * d + kernel["phi"]) / (1.0 + kernel["beta"] * (d**kernel["eta"]))
    w_raw = (np.abs(kvec) ** params["p_amp"]) * (d**params["r_dist"])
    w = w_raw / np.sum(w_raw)
    base_score = feats @ w

    lag_s = dff["best_lag_ms"].to_numpy(dtype=float) * 1e-3
    lag_phase_sin = np.sin(kernel["omega"] * lag_s + kernel["phi"])
    lag_phase_cos = np.cos(kernel["omega"] * lag_s + kernel["phi"])
    corr0 = dff["corr_at_0ms"].to_numpy(dtype=float)
    corr10 = dff["corr_at_10ms"].to_numpy(dtype=float)
    mean_abs = dff["mean_abs_corr"].to_numpy(dtype=float)

    pair_sign = np.where(pairs == 1, 1.0, np.where(pairs == 2, -1.0, 0.0))
    c1_raw = pair_sign * lag_phase_sin * (corr0 - corr10)
    c3_raw = pair_sign * lag_phase_cos * (mean_abs + np.abs(corr0) + np.abs(corr10))
    c4_raw = pair_sign * (lag_phase_sin * lag_phase_cos) * (np.abs(corr0) + np.abs(corr10) + mean_abs)

    c1 = c1_raw / max(float(np.std(c1_raw)), 1e-12)
    c3 = c3_raw / max(float(np.std(c3_raw)), 1e-12)
    c4 = c4_raw / max(float(np.std(c4_raw)), 1e-12)
    return base_score, c1, c3, c4, pairs


def compress_score(s_raw: np.ndarray, gamma_c: float) -> np.ndarray:
    med = float(np.median(s_raw))
    std = max(float(np.std(s_raw)), 1e-12)
    z = (s_raw - med) / std
    zc = np.sign(z) * (np.abs(z) ** gamma_c)
    return med + std * zc


def eval_fold(
    dff: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
    thr: Dict[str, float],
    xi1: float,
    xi3: float,
    xi4: float,
    gamma_c: float,
    real_iid_boot: int,
    null_trials: int,
    null_boot: int,
    seed: int,
) -> Dict[str, float]:
    base, c1, c3, c4, pairs = build_components(dff, kernel, params)
    s_raw = base + xi1 * c1 + xi3 * c3 + xi4 * c4
    s = compress_score(s_raw, gamma_c)
    s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]

    real_iid = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=real_iid_boot, seed=seed + 1)

    ctrl_idx = np.where(pairs != 0)[0]
    n_ctrl = len(ctrl_idx)
    n_plus = int(np.sum(pairs == 1))
    rng = np.random.default_rng(seed + 100)
    rates = []
    for t in range(null_trials):
        signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
        rng.shuffle(signs)
        rand_sign = np.zeros_like(pairs, dtype=float)
        rand_sign[ctrl_idx] = signs

        c1n_raw = rand_sign * c1
        c3n_raw = rand_sign * c3
        c4n_raw = rand_sign * c4
        c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
        c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
        c4n = c4n_raw / max(float(np.std(c4n_raw)), 1e-12)
        s_n_raw = base + xi1 * c1n + xi3 * c3n + xi4 * c4n
        s_n = compress_score(s_n_raw, gamma_c)
        s_hl_n, s_hv_n, s_lv_n = s_n[pairs == 0], s_n[pairs == 1], s_n[pairs == 2]
        rb = bootstrap_pass_rate(s_hl_n, s_hv_n, s_lv_n, thr, n_boot=null_boot, seed=seed + 200 + t)
        rates.append(float(rb))

    arr = np.array(rates, dtype=float)
    return {
        "real_iid": float(real_iid),
        "null_random_mean": float(np.mean(arr)),
        "null_random_p90": float(np.quantile(arr, 0.90)),
    }


def constrained_adv_rate_compressed(
    dff: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
    thr: Dict[str, float],
    xi1: float,
    xi3: float,
    xi4: float,
    gamma_c: float,
    max_flips: int,
    seed: int,
) -> float:
    base, c1, c3, c4, pairs = build_components(dff, kernel, params)
    ctrl_idx = np.where(pairs != 0)[0]
    n_plus = int(np.sum(pairs == 1))
    order = np.argsort(dff.iloc[ctrl_idx]["window_idx"].to_numpy(dtype=int))
    ordered_ctrl = ctrl_idx[order]
    t = xi1 * c1 + xi3 * c3 + xi4 * c4
    signs_ord = optimal_sign_sequence_with_flip_budget(t[ordered_ctrl], n_plus=n_plus, max_flips=max_flips)

    rand_sign = np.zeros_like(pairs, dtype=float)
    rand_sign[ordered_ctrl] = signs_ord
    c1n_raw = rand_sign * c1
    c3n_raw = rand_sign * c3
    c4n_raw = rand_sign * c4
    c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
    c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
    c4n = c4n_raw / max(float(np.std(c4n_raw)), 1e-12)
    s_n_raw = base + xi1 * c1n + xi3 * c3n + xi4 * c4n
    s_n = compress_score(s_n_raw, gamma_c)
    s_hl_n, s_hv_n, s_lv_n = s_n[pairs == 0], s_n[pairs == 1], s_n[pairs == 2]
    return float(bootstrap_pass_rate(s_hl_n, s_hv_n, s_lv_n, thr, n_boot=FULL_ADV_BOOT, seed=seed))


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    r1992 = json.loads(IN_QW1992.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]
    center = r1992["best"]

    rows = []
    d13 = np.array([-0.00007, -0.000035, 0.0, 0.000035, 0.00007], dtype=float)
    d4 = np.array([-0.00012, -0.00006, 0.0, 0.00006, 0.00012], dtype=float)
    gammas = np.array([0.82, 0.88, 0.94, 0.98, 1.0], dtype=float)
    for a in d13:
        for b in d13:
            for d in d4:
                for g in gammas:
                    rows.append(
                        {
                            "xi1": float(np.clip(float(center["xi1"]) + a, 0.0001, 0.0032)),
                            "xi3": float(np.clip(float(center["xi3"]) + b, -0.0002, 0.0032)),
                            "xi4": float(np.clip(float(center["xi4"]) + d, -0.0032, 0.0032)),
                            "gamma_c": float(g),
                        }
                    )
    rng = np.random.default_rng(230000)
    for _ in range(120):
        rows.append(
            {
                "xi1": float(np.clip(rng.normal(float(center["xi1"]), 0.00006), 0.0001, 0.0032)),
                "xi3": float(np.clip(rng.normal(float(center["xi3"]), 0.00006), -0.0002, 0.0032)),
                "xi4": float(np.clip(rng.normal(float(center["xi4"]), 0.00009), -0.0032, 0.0032)),
                "gamma_c": float(np.clip(rng.normal(0.94, 0.06), 0.78, 1.02)),
            }
        )

    seen = set()
    candidates = []
    for r in rows:
        k = (
            round(float(r["xi1"]), 12),
            round(float(r["xi3"]), 12),
            round(float(r["xi4"]), 12),
            round(float(r["gamma_c"]), 6),
        )
        if k in seen:
            continue
        seen.add(k)
        candidates.append({"xi1": k[0], "xi3": k[1], "xi4": k[2], "gamma_c": k[3]})

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    fast_rows = []
    total = len(candidates)
    for i, c in enumerate(candidates):
        xi1, xi3, xi4, gc = float(c["xi1"]), float(c["xi3"]), float(c["xi4"]), float(c["gamma_c"])
        fr = []
        for f, dff in enumerate(fold_dfs):
            row = eval_fold(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                gamma_c=gc,
                real_iid_boot=FAST_REAL_IID_BOOT,
                null_trials=FAST_NULL_TRIALS,
                null_boot=FAST_NULL_BOOT,
                seed=231000 + i * 10 + f,
            )
            fr.append({"fold": f, **row})

        min_real = float(min(r["real_iid"] for r in fr))
        max_null = float(max(r["null_random_p90"] for r in fr))
        hard_margin = float(min(min_real - TARGET_REAL_IID_MIN, TARGET_NULL_RANDOM_P90_MAX - max_null))
        aux = float(np.mean([r["real_iid"] for r in fr]) - np.mean([r["null_random_p90"] for r in fr]))
        fast_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "xi4": xi4,
                "gamma_c": gc,
                "fold_fast": fr,
                "min_real_iid_fast": min_real,
                "max_null_random_p90_fast": max_null,
                "hard_margin_fast": hard_margin,
                "aux_fast": aux,
            }
        )
        if (i + 1) % 60 == 0:
            print(f"[QW-1993] fast search progress: {i + 1}/{total}", flush=True)

    fast_rows.sort(key=lambda x: (x["hard_margin_fast"], x["aux_fast"]), reverse=True)
    shortlist = fast_rows[:SHORTLIST_SIZE]
    print(f"[QW-1993] shortlist size: {len(shortlist)}", flush=True)

    final_rows = []
    for i, c in enumerate(shortlist):
        xi1, xi3, xi4, gc = float(c["xi1"]), float(c["xi3"]), float(c["xi4"]), float(c["gamma_c"])
        fr = []
        for f, dff in enumerate(fold_dfs):
            base = eval_fold(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                gamma_c=gc,
                real_iid_boot=FULL_REAL_IID_BOOT,
                null_trials=FULL_NULL_TRIALS,
                null_boot=FULL_NULL_BOOT,
                seed=235000 + i * 10 + f,
            )
            comp = build_components(dff, kernel, params)
            base_score, c1, c3, c4, pairs = comp
            s_raw = base_score + xi1 * c1 + xi3 * c3 + xi4 * c4
            s = compress_score(s_raw, gc)
            s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]
            real_block = bootstrap_pass_rate_block(
                s_hl=s_hl,
                s_hv=s_hv,
                s_lv=s_lv,
                thr=thr,
                n_boot=FULL_REAL_BLOCK_BOOT,
                seed=236000 + i * 10 + f,
                block_len=FULL_REAL_BLOCK_LEN,
            )
            adv = constrained_adv_rate_compressed(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                gamma_c=gc,
                max_flips=FULL_ADV_FLIP_BUDGET,
                seed=237000 + i * 10 + f,
            )
            fr.append(
                {
                    "fold": f,
                    "real_iid_full": base["real_iid"],
                    "real_block_full": real_block,
                    "null_random_mean_full": base["null_random_mean"],
                    "null_random_p90_full": base["null_random_p90"],
                    "adv_constrained_full": adv,
                }
            )

        min_real_iid = float(min(r["real_iid_full"] for r in fr))
        min_real_block = float(min(r["real_block_full"] for r in fr))
        max_null = float(max(r["null_random_p90_full"] for r in fr))
        max_adv = float(max(r["adv_constrained_full"] for r in fr))
        hard_margin = float(min(min_real_iid - TARGET_REAL_IID_MIN, TARGET_NULL_RANDOM_P90_MAX - max_null))
        hard_pass = bool(
            min_real_iid >= TARGET_REAL_IID_MIN
            and min_real_block >= TARGET_REAL_BLOCK_MIN
            and max_null <= TARGET_NULL_RANDOM_P90_MAX
            and max_adv <= 0.45
        )
        aux = float(np.mean([r["real_iid_full"] for r in fr]) - np.mean([r["null_random_p90_full"] for r in fr]))
        final_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "xi4": xi4,
                "gamma_c": gc,
                "hard_pass": hard_pass,
                "min_real_iid_full": min_real_iid,
                "min_real_block_full": min_real_block,
                "max_null_random_p90_full": max_null,
                "max_adv_constrained_full": max_adv,
                "hard_margin_full": hard_margin,
                "aux_full": aux,
                "fold_results": fr,
            }
        )

    final_rows.sort(
        key=lambda x: (
            int(x["hard_pass"]),
            x["hard_margin_full"],
            x["aux_full"],
        ),
        reverse=True,
    )
    best = final_rows[0]
    verdict = "COMPRESSION_HARD_PASS" if best["hard_pass"] else "COMPRESSION_HARD_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1992.name],
        "search_config": {
            "total_candidates": total,
            "shortlist_size": SHORTLIST_SIZE,
            "target_real_iid_min": TARGET_REAL_IID_MIN,
            "target_null_random_p90_max": TARGET_NULL_RANDOM_P90_MAX,
            "target_real_block_min": TARGET_REAL_BLOCK_MIN,
            "fast_real_iid_boot": FAST_REAL_IID_BOOT,
            "fast_null_trials": FAST_NULL_TRIALS,
            "fast_null_boot": FAST_NULL_BOOT,
            "full_real_iid_boot": FULL_REAL_IID_BOOT,
            "full_real_block_boot": FULL_REAL_BLOCK_BOOT,
            "full_null_trials": FULL_NULL_TRIALS,
            "full_null_boot": FULL_NULL_BOOT,
            "full_adv_flip_budget": FULL_ADV_FLIP_BUDGET,
            "full_adv_boot": FULL_ADV_BOOT,
        },
        "best": best,
        "top8": final_rows,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1993: GLOBAL SCORE COMPRESSION HARD-LOCK SEARCH",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Best Candidate",
        f"- xi1/xi3/xi4/gamma_c: {best['xi1']:.6f}/{best['xi3']:.6f}/{best['xi4']:.6f}/{best['gamma_c']:.3f}",
        f"- hard_pass: {best['hard_pass']}",
        f"- min_real_iid_full: {100.0 * best['min_real_iid_full']:.2f}%",
        f"- min_real_block_full: {100.0 * best['min_real_block_full']:.2f}%",
        f"- max_null_random_p90_full: {100.0 * best['max_null_random_p90_full']:.2f}%",
        f"- max_adv_constrained_full: {100.0 * best['max_adv_constrained_full']:.2f}%",
        f"- hard_margin_full: {100.0 * best['hard_margin_full']:.2f} pp",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1993] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1993] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1993] verdict={verdict}")


if __name__ == "__main__":
    main()

