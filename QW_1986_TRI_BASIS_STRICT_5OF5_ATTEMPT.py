#!/usr/bin/env python3
"""
QW-1986: Tri-basis frozen-operator strict 5/5 attempt.
Adds one mixed-phase coupling channel (xi4) to reduce worst-fold null leakage.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import ROOT, bootstrap_pass_rate, gw_flags, gw_metrics


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1985 = ROOT / "report_qw1985_strict_5of5_local_push.json"
OUT_JSON = ROOT / "report_qw1986_tri_basis_strict_5of5_attempt.json"
OUT_MD = ROOT / "RAPORT_QW1986_TRI_BASIS_STRICT_5OF5_ATTEMPT.md"

REAL_PASS_MIN = 0.85
NULL_P90_PASS_MAX = 0.40

FAST_REAL_BOOT = 90
FAST_NULL_TRIALS = 6
FAST_NULL_BOOT = 30

FULL_REAL_BOOT = 1000
FULL_NULL_TRIALS = 14
FULL_NULL_BOOT = 70
SHORTLIST_SIZE = 14


def build_fold_channels_tri(
    df_fold: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
    xi1: float,
    xi3: float,
    xi4: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
    # New mixed-phase odd channel (still signed by control-pair orientation).
    c4_raw = pair_sign * (lag_phase_sin * lag_phase_cos) * (np.abs(corr0) + np.abs(corr10) + mean_abs)

    c1 = c1_raw / max(float(np.std(c1_raw)), 1e-12)
    c3 = c3_raw / max(float(np.std(c3_raw)), 1e-12)
    c4 = c4_raw / max(float(np.std(c4_raw)), 1e-12)

    s = base_score + xi1 * c1 + xi3 * c3 + xi4 * c4
    return s[pairs == 0], s[pairs == 1], s[pairs == 2], pairs, c1, c3, c4


def fold_null_stats_tri(
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
    pair_sign = np.where(pairs == 1, 1.0, np.where(pairs == 2, -1.0, 0.0))

    random_rates: List[float] = []
    for i in range(n_trials):
        signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
        rng.shuffle(signs)
        rand_sign = np.zeros_like(pair_sign)
        rand_sign[ctrl_idx] = signs

        c1n_raw = rand_sign * c1
        c3n_raw = rand_sign * c3
        c4n_raw = rand_sign * c4
        c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
        c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
        c4n = c4n_raw / max(float(np.std(c4n_raw)), 1e-12)
        s = s_base + xi1 * c1n + xi3 * c3n + xi4 * c4n
        s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]
        rb = bootstrap_pass_rate(
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            thr=thr,
            n_boot=n_boot,
            seed=seed + 300 + i,
        )
        random_rates.append(float(rb))

    arr = np.array(random_rates, dtype=float)
    return float(np.mean(arr)), float(np.quantile(arr, 0.90))


def eval_on_fold(
    dff: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
    thr: Dict[str, float],
    xi1: float,
    xi3: float,
    xi4: float,
    real_boot: int,
    null_trials: int,
    null_boot: int,
    seed: int,
) -> Dict[str, float]:
    s_hl, s_hv, s_lv, pairs, c1, c3, c4 = build_fold_channels_tri(dff, kernel, params, xi1, xi3, xi4)
    g = gw_metrics(s_hl, s_hv, s_lv)
    flags = gw_flags(g, thr)
    real_rate = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=real_boot, seed=seed + 1)

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
        seed=seed + 70,
        n_trials=null_trials,
        n_boot=null_boot,
    )
    return {
        "real_rate": float(real_rate),
        "null_mean": float(null_mean),
        "null_p90": float(null_p90),
        "det_flags": flags,
    }


def unique_rows(rows: List[Dict[str, float]]) -> List[Dict[str, float]]:
    seen = set()
    out = []
    for r in rows:
        k = (round(float(r["xi1"]), 12), round(float(r["xi3"]), 12), round(float(r["xi4"]), 12))
        if k in seen:
            continue
        seen.add(k)
        out.append({"xi1": float(k[0]), "xi3": float(k[1]), "xi4": float(k[2])})
    return out


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    r1985 = json.loads(IN_QW1985.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    xi1_c = float(r1985["best"]["xi1"])
    xi3_c = float(r1985["best"]["xi3"])
    xi4_c = 0.0

    d13 = np.array([-0.00012, -0.00006, 0.0, 0.00006, 0.00012], dtype=float)
    d4 = np.array([-0.0008, -0.0004, 0.0, 0.0004, 0.0008], dtype=float)

    rows = []
    for a in d13:
        for b in d13:
            for c in d4:
                rows.append(
                    {
                        "xi1": float(np.clip(xi1_c + a, 0.0001, 0.0028)),
                        "xi3": float(np.clip(xi3_c + b, -0.0002, 0.0028)),
                        "xi4": float(np.clip(xi4_c + c, -0.0025, 0.0025)),
                    }
                )
    rng = np.random.default_rng(160000)
    for _ in range(120):
        rows.append(
            {
                "xi1": float(np.clip(rng.normal(xi1_c, 0.00010), 0.0001, 0.0028)),
                "xi3": float(np.clip(rng.normal(xi3_c, 0.00010), -0.0002, 0.0028)),
                "xi4": float(np.clip(rng.normal(0.0, 0.0006), -0.0025, 0.0025)),
            }
        )
    candidates = unique_rows(rows)

    fast_rows = []
    total = len(candidates)
    for i, c in enumerate(candidates):
        xi1 = float(c["xi1"])
        xi3 = float(c["xi3"])
        xi4 = float(c["xi4"])
        fold_fast = []
        for f, dff in enumerate(fold_dfs):
            row = eval_on_fold(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                real_boot=FAST_REAL_BOOT,
                null_trials=FAST_NULL_TRIALS,
                null_boot=FAST_NULL_BOOT,
                seed=161000 + i * 10 + f,
            )
            fold_fast.append(
                {
                    "fold": f,
                    "real_fast": row["real_rate"],
                    "null_p90_fast": row["null_p90"],
                    "null_mean_fast": row["null_mean"],
                    "det_flags": row["det_flags"],
                }
            )
        min_real = float(min(r["real_fast"] for r in fold_fast))
        max_null = float(max(r["null_p90_fast"] for r in fold_fast))
        strict_margin = float(min(min_real - REAL_PASS_MIN, NULL_P90_PASS_MAX - max_null))
        strict_fast = bool(strict_margin >= 0.0)
        aux = float(np.mean([r["real_fast"] for r in fold_fast]) - np.mean([r["null_p90_fast"] for r in fold_fast]))
        fast_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "xi4": xi4,
                "fold_fast": fold_fast,
                "min_real_fast": min_real,
                "max_null_p90_fast": max_null,
                "strict_margin_fast": strict_margin,
                "strict_fast": strict_fast,
                "aux_score_fast": aux,
            }
        )
        if (i + 1) % 40 == 0:
            print(f"[QW-1986] fast search progress: {i + 1}/{total}", flush=True)

    fast_rows.sort(
        key=lambda x: (int(x["strict_fast"]), x["strict_margin_fast"], x["aux_score_fast"]),
        reverse=True,
    )
    shortlist = fast_rows[:SHORTLIST_SIZE]
    print(f"[QW-1986] shortlist size: {len(shortlist)}", flush=True)

    final_rows = []
    for i, c in enumerate(shortlist):
        xi1 = float(c["xi1"])
        xi3 = float(c["xi3"])
        xi4 = float(c["xi4"])
        fold_full = []
        for f, dff in enumerate(fold_dfs):
            row = eval_on_fold(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                real_boot=FULL_REAL_BOOT,
                null_trials=FULL_NULL_TRIALS,
                null_boot=FULL_NULL_BOOT,
                seed=165000 + i * 10 + f,
            )
            fold_pass = bool(row["real_rate"] >= REAL_PASS_MIN and row["null_p90"] <= NULL_P90_PASS_MAX)
            fold_full.append(
                {
                    "fold": f,
                    "real_full": row["real_rate"],
                    "null_mean_full": row["null_mean"],
                    "null_p90_full": row["null_p90"],
                    "det_flags": row["det_flags"],
                    "fold_pass": fold_pass,
                }
            )
        pass_count = int(sum(int(r["fold_pass"]) for r in fold_full))
        min_real = float(min(r["real_full"] for r in fold_full))
        max_null = float(max(r["null_p90_full"] for r in fold_full))
        strict_margin = float(min(min_real - REAL_PASS_MIN, NULL_P90_PASS_MAX - max_null))
        aux = float(np.mean([r["real_full"] for r in fold_full]) - np.mean([r["null_p90_full"] for r in fold_full]))
        final_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "xi4": xi4,
                "pass_count": pass_count,
                "min_real_full": min_real,
                "max_null_p90_full": max_null,
                "strict_margin_full": strict_margin,
                "aux_score_full": aux,
                "fold_results": fold_full,
            }
        )

    final_rows.sort(
        key=lambda x: (x["pass_count"], x["strict_margin_full"], x["aux_score_full"]),
        reverse=True,
    )
    best = final_rows[0]

    verdict = (
        "TRI_BASIS_STRICT_5OF5_PASS"
        if best["pass_count"] == 5
        else "TRI_BASIS_STRICT_NEAR_MISS"
        if best["pass_count"] == 4
        else "TRI_BASIS_STRICT_FAIL"
    )
    required = (
        "LOCK_TRI_BASIS_OPERATOR_AND_RUN_HARD_TEMPORAL_EXTERNAL_TESTS"
        if best["pass_count"] == 5
        else "BORDERLINE_PERSISTS_REQUIRES_STRONGER_PHYSICAL_CONSTRAINT_OR_NEW_INVARIANT"
        if best["pass_count"] == 4
        else "TRI_BASIS_INSUFFICIENT_EXTEND_MODEL_CLASS"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1985.name],
        "search_config": {
            "grid_points": len(d13) * len(d13) * len(d4),
            "random_points": 120,
            "total_candidates": total,
            "shortlist_size": SHORTLIST_SIZE,
            "real_pass_min": REAL_PASS_MIN,
            "null_p90_pass_max": NULL_P90_PASS_MAX,
            "fast_real_boot": FAST_REAL_BOOT,
            "fast_null_trials": FAST_NULL_TRIALS,
            "fast_null_boot": FAST_NULL_BOOT,
            "full_real_boot": FULL_REAL_BOOT,
            "full_null_trials": FULL_NULL_TRIALS,
            "full_null_boot": FULL_NULL_BOOT,
        },
        "center": {"xi1": xi1_c, "xi3": xi3_c, "xi4": xi4_c},
        "best": best,
        "top14": final_rows,
        "verdict": verdict,
        "required_next_step": required,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1986: TRI-BASIS STRICT 5/5 ATTEMPT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Best Candidate",
        f"- xi1/xi3/xi4: {best['xi1']:.6f}/{best['xi3']:.6f}/{best['xi4']:.6f}",
        f"- pass_count: {best['pass_count']}/5",
        f"- min_real_full: {100.0 * best['min_real_full']:.2f}%",
        f"- max_null_p90_full: {100.0 * best['max_null_p90_full']:.2f}%",
        f"- strict_margin_full: {100.0 * best['strict_margin_full']:.2f} pp",
        "",
        "## Required Next Step",
        f"- {required}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1986] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1986] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1986] verdict={verdict} pass_count={best['pass_count']}/5")


if __name__ == "__main__":
    main()

