#!/usr/bin/env python3
"""
QW-1984: Min-max fold refinement for a single frozen operator.
Goal: improve worst-fold null leakage without fold-specific retune.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import (
    ROOT,
    bootstrap_pass_rate,
    build_fold_channels,
    gw_flags,
    gw_metrics,
)


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1983 = ROOT / "report_qw1983_fold_robust_operator_search.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
OUT_JSON = ROOT / "report_qw1984_fold_minmax_nullp90_refinement.json"
OUT_MD = ROOT / "RAPORT_QW1984_FOLD_MINMAX_NULLP90_REFINEMENT.md"

N_SEEDS = 4
N_LOCAL_PER_SEED = 70
SIGMA_XI1 = 0.00012
SIGMA_XI3 = 0.00012

FAST_REAL_BOOT = 90
FAST_NULL_TRIALS = 6
FAST_NULL_BOOT = 30

FULL_REAL_BOOT = 800
FULL_NULL_TRIALS = 12
FULL_NULL_BOOT = 60
SHORTLIST_SIZE = 10

REAL_PASS_MIN = 0.85
NULL_P90_PASS_MAX = 0.40


def fold_null_stats(
    s_base: np.ndarray,
    pairs: np.ndarray,
    c1: np.ndarray,
    c3: np.ndarray,
    xi1: float,
    xi3: float,
    thr: Dict[str, float],
    seed: int,
    n_trials: int,
    n_boot: int,
) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    ctrl_idx = np.where(pairs != 0)[0]
    n_ctrl = len(ctrl_idx)
    n_plus = int(np.sum(pairs == 1))
    pair_sign = np.where(pairs == 1, 1.0, np.where(pairs == 2, -1.0, 0.0))

    random_rates: List[float] = []
    for i in range(n_trials):
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
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            thr=thr,
            n_boot=n_boot,
            seed=seed + 200 + i,
        )
        random_rates.append(float(rb))

    arr = np.array(random_rates, dtype=float)
    return float(np.mean(arr)), float(np.quantile(arr, 0.90))


def eval_candidate_on_fold(
    dff: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
    thr: Dict[str, float],
    xi1: float,
    xi3: float,
    real_boot: int,
    null_trials: int,
    null_boot: int,
    seed: int,
) -> Dict[str, float]:
    s_hl, s_hv, s_lv, pairs, c1, c3 = build_fold_channels(dff, kernel, params, xi1, xi3)
    g = gw_metrics(s_hl, s_hv, s_lv)
    flags = gw_flags(g, thr)
    real_rate = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=real_boot, seed=seed + 1)

    pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
    pairs_vec = dff["pair"].map(pair_map).to_numpy(dtype=int)
    s_full = np.zeros(len(pairs_vec), dtype=float)
    s_full[pairs_vec == 0] = s_hl
    s_full[pairs_vec == 1] = s_hv
    s_full[pairs_vec == 2] = s_lv
    s_base = s_full - xi1 * c1 - xi3 * c3

    null_mean, null_p90 = fold_null_stats(
        s_base=s_base,
        pairs=pairs_vec,
        c1=c1,
        c3=c3,
        xi1=xi1,
        xi3=xi3,
        thr=thr,
        seed=seed + 50,
        n_trials=null_trials,
        n_boot=null_boot,
    )

    return {
        "real_rate": float(real_rate),
        "null_mean": float(null_mean),
        "null_p90": float(null_p90),
        "det_flags": flags,
    }


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1983 = json.loads(IN_QW1983.read_text(encoding="utf-8"))
    thr = json.loads(IN_QW1969.read_text(encoding="utf-8"))["thresholds"]
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    seed_rows = [r for r in r1983["top8"] if int(r["pass_count"]) >= 2][:N_SEEDS]
    if len(seed_rows) < N_SEEDS:
        seed_rows = r1983["top8"][:N_SEEDS]

    rng = np.random.default_rng(140000)
    candidate_pool: List[Dict[str, float]] = []
    for s_idx, seed in enumerate(seed_rows):
        xi1_0 = float(seed["xi1"])
        xi3_0 = float(seed["xi3"])
        candidate_pool.append({"xi1": xi1_0, "xi3": xi3_0, "seed_idx": s_idx, "is_seed": True})
        for _ in range(N_LOCAL_PER_SEED):
            candidate_pool.append(
                {
                    "xi1": float(np.clip(rng.normal(xi1_0, SIGMA_XI1), 0.0001, 0.0026)),
                    "xi3": float(np.clip(rng.normal(xi3_0, SIGMA_XI3), -0.0002, 0.0024)),
                    "seed_idx": s_idx,
                    "is_seed": False,
                }
            )

    fast_rows = []
    total = len(candidate_pool)
    for i, c in enumerate(candidate_pool):
        xi1 = float(c["xi1"])
        xi3 = float(c["xi3"])
        fold_rows = []
        for f, dff in enumerate(fold_dfs):
            row = eval_candidate_on_fold(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                real_boot=FAST_REAL_BOOT,
                null_trials=FAST_NULL_TRIALS,
                null_boot=FAST_NULL_BOOT,
                seed=141000 + i * 10 + f,
            )
            fold_rows.append(
                {
                    "fold": f,
                    "real_fast": row["real_rate"],
                    "null_p90_fast": row["null_p90"],
                    "null_mean_fast": row["null_mean"],
                    "det_flags": row["det_flags"],
                }
            )

        min_real = float(min(r["real_fast"] for r in fold_rows))
        max_null_p90 = float(max(r["null_p90_fast"] for r in fold_rows))
        mean_real = float(np.mean([r["real_fast"] for r in fold_rows]))
        mean_null = float(np.mean([r["null_p90_fast"] for r in fold_rows]))
        strict_fast = bool(min_real >= REAL_PASS_MIN and max_null_p90 <= NULL_P90_PASS_MAX)
        # Min-max score over hardest fold.
        minmax_margin = float(min(min_real - REAL_PASS_MIN, NULL_P90_PASS_MAX - max_null_p90))
        aux_score = float(mean_real - mean_null)
        fast_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "seed_idx": int(c["seed_idx"]),
                "is_seed": bool(c["is_seed"]),
                "fold_fast": fold_rows,
                "min_real_fast": min_real,
                "max_null_p90_fast": max_null_p90,
                "mean_real_fast": mean_real,
                "mean_null_p90_fast": mean_null,
                "strict_fast": strict_fast,
                "minmax_margin_fast": minmax_margin,
                "aux_score_fast": aux_score,
            }
        )
        if (i + 1) % 40 == 0:
            print(f"[QW-1984] fast search progress: {i + 1}/{total}", flush=True)

    fast_rows.sort(
        key=lambda x: (
            int(x["strict_fast"]),
            x["minmax_margin_fast"],
            x["aux_score_fast"],
        ),
        reverse=True,
    )
    shortlist = fast_rows[:SHORTLIST_SIZE]
    print(f"[QW-1984] shortlist size: {len(shortlist)}", flush=True)

    final_rows = []
    for i, cand in enumerate(shortlist):
        xi1 = float(cand["xi1"])
        xi3 = float(cand["xi3"])
        folds_full = []
        for f, dff in enumerate(fold_dfs):
            row = eval_candidate_on_fold(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                real_boot=FULL_REAL_BOOT,
                null_trials=FULL_NULL_TRIALS,
                null_boot=FULL_NULL_BOOT,
                seed=145000 + i * 10 + f,
            )
            fold_pass = bool(row["real_rate"] >= REAL_PASS_MIN and row["null_p90"] <= NULL_P90_PASS_MAX)
            folds_full.append(
                {
                    "fold": f,
                    "real_full": row["real_rate"],
                    "null_mean_full": row["null_mean"],
                    "null_p90_full": row["null_p90"],
                    "det_flags": row["det_flags"],
                    "fold_pass": fold_pass,
                }
            )

        pass_count = int(sum(int(r["fold_pass"]) for r in folds_full))
        min_real = float(min(r["real_full"] for r in folds_full))
        max_null_p90 = float(max(r["null_p90_full"] for r in folds_full))
        minmax_margin = float(min(min_real - REAL_PASS_MIN, NULL_P90_PASS_MAX - max_null_p90))
        aux_score = float(
            np.mean([r["real_full"] for r in folds_full]) - np.mean([r["null_p90_full"] for r in folds_full])
        )
        final_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "pass_count": pass_count,
                "min_real_full": min_real,
                "max_null_p90_full": max_null_p90,
                "minmax_margin_full": minmax_margin,
                "aux_score_full": aux_score,
                "fold_results": folds_full,
            }
        )

    final_rows.sort(
        key=lambda x: (
            x["pass_count"],
            x["minmax_margin_full"],
            x["aux_score_full"],
        ),
        reverse=True,
    )
    best = final_rows[0]

    verdict = (
        "FOLD_MINMAX_REFINEMENT_PASS_STRONG"
        if best["pass_count"] >= 4
        else "FOLD_MINMAX_REFINEMENT_PASS_PARTIAL"
        if best["pass_count"] >= 3
        else "FOLD_MINMAX_REFINEMENT_FAIL"
    )
    required = (
        "PROMOTE_TO_HARDER_TEMPORAL_AND_EXTERNAL_CONFIRMATORY_PREP"
        if best["pass_count"] >= 4
        else "REPEAT_REFINEMENT_OR_EXTEND_OPERATOR_BASIS"
        if best["pass_count"] >= 3
        else "EXTEND_OPERATOR_BASIS_AND_RESTART_SEARCH"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1983.name, IN_QW1969.name],
        "search_config": {
            "n_seeds": N_SEEDS,
            "n_local_per_seed": N_LOCAL_PER_SEED,
            "total_candidates": total,
            "shortlist_size": SHORTLIST_SIZE,
            "fast_real_boot": FAST_REAL_BOOT,
            "fast_null_trials": FAST_NULL_TRIALS,
            "fast_null_boot": FAST_NULL_BOOT,
            "full_real_boot": FULL_REAL_BOOT,
            "full_null_trials": FULL_NULL_TRIALS,
            "full_null_boot": FULL_NULL_BOOT,
            "real_pass_min": REAL_PASS_MIN,
            "null_p90_pass_max": NULL_P90_PASS_MAX,
        },
        "best": best,
        "top10": final_rows,
        "verdict": verdict,
        "required_next_step": required,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1984: FOLD MIN-MAX NULL-P90 REFINEMENT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Best Candidate",
        f"- xi1/xi3: {best['xi1']:.6f}/{best['xi3']:.6f}",
        f"- pass_count: {best['pass_count']}/5",
        f"- min_real_full: {100.0 * best['min_real_full']:.2f}%",
        f"- max_null_p90_full: {100.0 * best['max_null_p90_full']:.2f}%",
        f"- minmax_margin_full: {100.0 * best['minmax_margin_full']:.2f} pp",
        "",
        "## Required Next Step",
        f"- {required}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1984] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1984] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1984] verdict={verdict} pass_count={best['pass_count']}/5")


if __name__ == "__main__":
    main()

