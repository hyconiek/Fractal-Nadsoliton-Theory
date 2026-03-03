#!/usr/bin/env python3
"""
QW-1988: Hard lock stress-test for the QW-1987 tri-basis winner.
Checks IID and block stability plus random/adversarial null leakage per fold.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import ROOT, bootstrap_pass_rate
from QW_1986_TRI_BASIS_STRICT_5OF5_ATTEMPT import build_fold_channels_tri, fold_null_stats_tri


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1987 = ROOT / "report_qw1987_tri_basis_fold2_targeted_strict_push.json"
OUT_JSON = ROOT / "report_qw1988_tri_basis_hard_lock_stress.json"
OUT_MD = ROOT / "RAPORT_QW1988_TRI_BASIS_HARD_LOCK_STRESS.md"

REAL_IID_BOOT = 3000
REAL_BLOCK_BOOT = 1800
REAL_BLOCK_LEN = 10

NULL_RANDOM_TRIALS = 32
NULL_RANDOM_BOOT = 120
NULL_ADV_TRIALS = 24
NULL_ADV_BOOT = 120


def bootstrap_pass_rate_block(
    s_hl: np.ndarray,
    s_hv: np.ndarray,
    s_lv: np.ndarray,
    thr: Dict[str, float],
    n_boot: int,
    seed: int,
    block_len: int,
) -> float:
    rng = np.random.default_rng(seed)

    def resample_blocks(x: np.ndarray) -> np.ndarray:
        n = len(x)
        n_blocks = int(np.ceil(n / block_len))
        starts = rng.integers(0, max(1, n - block_len + 1), size=n_blocks, endpoint=False)
        chunks = [x[s : s + block_len] for s in starts]
        y = np.concatenate(chunks)
        return y[:n]

    pass_count = 0
    for _ in range(n_boot):
        b_hl = resample_blocks(s_hl)
        b_hv = resample_blocks(s_hv)
        b_lv = resample_blocks(s_lv)
        if _all_flags(b_hl, b_hv, b_lv, thr):
            pass_count += 1
    return float(pass_count / n_boot)


def _all_flags(s_hl: np.ndarray, s_hv: np.ndarray, s_lv: np.ndarray, thr: Dict[str, float]) -> bool:
    s_ctrl = np.concatenate([s_hv, s_lv])
    q90 = float(np.quantile(s_ctrl, 0.90))
    adv = float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90))
    sep = float(np.median(s_hl) - np.median(s_ctrl))
    gap = float(abs(np.median(s_hv) - np.median(s_lv)))
    auc = _rank_auc_pos_gt_neg(s_hl, s_ctrl)
    return bool(
        sep >= thr["gw_sep_min"]
        and adv >= thr["gw_adv_min"]
        and auc >= thr["gw_auc_min"]
        and gap <= thr["gw_control_gap_max"]
    )


def _rank_auc_pos_gt_neg(pos: np.ndarray, neg: np.ndarray) -> float:
    y = np.concatenate([np.ones(len(pos), dtype=int), np.zeros(len(neg), dtype=int)])
    s = np.concatenate([pos, neg])
    n1 = len(pos)
    n0 = len(neg)
    order = np.argsort(s)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    rs = float(np.sum(ranks[y == 1]))
    return float((rs - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def adversarial_null_stats_tri(
    s_base: np.ndarray,
    pairs: np.ndarray,
    c1: np.ndarray,
    c3: np.ndarray,
    c4: np.ndarray,
    xi1: float,
    xi3: float,
    xi4: float,
    thr: Dict[str, float],
    seed: int,
    n_trials: int,
    n_boot: int,
) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    ctrl_idx = np.where(pairs != 0)[0]
    n_ctrl = len(ctrl_idx)
    n_plus = int(np.sum(pairs == 1))

    # adversarial template: maximize channel score on controls
    t = xi1 * c1 + xi3 * c3 + xi4 * c4
    order_ctrl = ctrl_idx[np.argsort(t[ctrl_idx])]

    rates = []
    for i in range(n_trials):
        signs = -np.ones(n_ctrl, dtype=float)
        # assign +1 to strongest controls; add small random swaps to avoid single deterministic pattern
        plus_idx = order_ctrl[-n_plus:].copy()
        if n_plus > 0:
            n_swaps = max(1, n_plus // 12)
            swap_out = rng.choice(np.arange(n_plus), size=n_swaps, replace=False)
            swap_in = rng.choice(np.arange(n_ctrl - n_plus), size=n_swaps, replace=False)
            plus_idx[swap_out] = order_ctrl[swap_in]
        mask_plus = np.isin(order_ctrl, plus_idx)
        signs[mask_plus] = 1.0

        rand_sign = np.zeros_like(pairs, dtype=float)
        rand_sign[ctrl_idx] = signs

        c1n_raw = rand_sign * c1
        c3n_raw = rand_sign * c3
        c4n_raw = rand_sign * c4
        c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
        c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
        c4n = c4n_raw / max(float(np.std(c4n_raw)), 1e-12)
        s = s_base + xi1 * c1n + xi3 * c3n + xi4 * c4n
        s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]
        rb = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=n_boot, seed=seed + 500 + i)
        rates.append(float(rb))

    arr = np.array(rates, dtype=float)
    return float(np.mean(arr)), float(np.quantile(arr, 0.90))


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    r1987 = json.loads(IN_QW1987.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]
    cand = r1987["best"]
    xi1 = float(cand["xi1"])
    xi3 = float(cand["xi3"])
    xi4 = float(cand["xi4"])

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    fold_rows = []
    for f, dff in enumerate(fold_dfs):
        s_hl, s_hv, s_lv, pairs, c1, c3, c4 = build_fold_channels_tri(dff, kernel, params, xi1, xi3, xi4)
        real_iid = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=REAL_IID_BOOT, seed=181000 + f)
        real_block = bootstrap_pass_rate_block(
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            thr=thr,
            n_boot=REAL_BLOCK_BOOT,
            seed=182000 + f,
            block_len=REAL_BLOCK_LEN,
        )

        pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
        pairs_vec = dff["pair"].map(pair_map).to_numpy(dtype=int)
        s_full = np.zeros(len(pairs_vec), dtype=float)
        s_full[pairs_vec == 0] = s_hl
        s_full[pairs_vec == 1] = s_hv
        s_full[pairs_vec == 2] = s_lv
        s_base = s_full - xi1 * c1 - xi3 * c3 - xi4 * c4

        null_mean, null_p90 = fold_null_stats_tri(
            s_base=s_base,
            pairs=pairs_vec,
            c1=c1,
            c3=c3,
            c4=c4,
            xi1=xi1,
            xi3=xi3,
            xi4=xi4,
            thr=thr,
            seed=183000 + f,
            n_trials=NULL_RANDOM_TRIALS,
            n_boot=NULL_RANDOM_BOOT,
        )
        adv_mean, adv_p90 = adversarial_null_stats_tri(
            s_base=s_base,
            pairs=pairs_vec,
            c1=c1,
            c3=c3,
            c4=c4,
            xi1=xi1,
            xi3=xi3,
            xi4=xi4,
            thr=thr,
            seed=184000 + f,
            n_trials=NULL_ADV_TRIALS,
            n_boot=NULL_ADV_BOOT,
        )

        fold_rows.append(
            {
                "fold": f,
                "real_iid": real_iid,
                "real_block": real_block,
                "null_random_mean": null_mean,
                "null_random_p90": null_p90,
                "null_adv_mean": adv_mean,
                "null_adv_p90": adv_p90,
            }
        )
        print(f"[QW-1988] fold {f} done", flush=True)

    real_iid_min = float(min(r["real_iid"] for r in fold_rows))
    real_block_min = float(min(r["real_block"] for r in fold_rows))
    null_random_p90_max = float(max(r["null_random_p90"] for r in fold_rows))
    null_adv_mean_max = float(max(r["null_adv_mean"] for r in fold_rows))
    null_adv_p90_max = float(max(r["null_adv_p90"] for r in fold_rows))

    hard_pass = bool(
        real_iid_min >= 0.95
        and real_block_min >= 0.90
        and null_random_p90_max <= 0.40
        and null_adv_mean_max <= 0.45
    )
    verdict = "HARD_LOCK_READY_INTERNAL" if hard_pass else "HARD_LOCK_NOT_READY"
    required = (
        "PROCEED_TO_EXTERNAL_CONFIRMATORY_SUITE"
        if hard_pass
        else "REWORK_OPERATOR_OR_NULL_GUARDS_BEFORE_EXTERNAL_STAGE"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1987.name],
        "candidate": {"xi1": xi1, "xi3": xi3, "xi4": xi4},
        "stress_config": {
            "real_iid_boot": REAL_IID_BOOT,
            "real_block_boot": REAL_BLOCK_BOOT,
            "real_block_len": REAL_BLOCK_LEN,
            "null_random_trials": NULL_RANDOM_TRIALS,
            "null_random_boot": NULL_RANDOM_BOOT,
            "null_adv_trials": NULL_ADV_TRIALS,
            "null_adv_boot": NULL_ADV_BOOT,
        },
        "fold_results": fold_rows,
        "aggregate": {
            "real_iid_min": real_iid_min,
            "real_block_min": real_block_min,
            "null_random_p90_max": null_random_p90_max,
            "null_adv_mean_max": null_adv_mean_max,
            "null_adv_p90_max": null_adv_p90_max,
        },
        "verdict": verdict,
        "required_next_step": required,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1988: TRI-BASIS HARD LOCK STRESS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Aggregate",
        f"- real_iid_min: {100.0 * real_iid_min:.2f}%",
        f"- real_block_min: {100.0 * real_block_min:.2f}%",
        f"- null_random_p90_max: {100.0 * null_random_p90_max:.2f}%",
        f"- null_adv_mean_max: {100.0 * null_adv_mean_max:.2f}%",
        f"- null_adv_p90_max: {100.0 * null_adv_p90_max:.2f}%",
        "",
        "## Required Next Step",
        f"- {required}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1988] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1988] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1988] verdict={verdict}")


if __name__ == "__main__":
    main()

