#!/usr/bin/env python3
"""
QW-1999: Bounded coupling fold-2 guarded robust hard-lock search.
Adds one global saturation scale kappa_t for the signed coupling term:
  t = xi1*c1 + xi3*c3 + xi4*c4,
  t_eff = clip(t, -kappa_t*std(t), +kappa_t*std(t)).
Readout remains globally monotonic (power compression gamma_c) and frozen-kernel.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import ROOT
from QW_1988_TRI_BASIS_HARD_LOCK_STRESS import bootstrap_pass_rate_block
from QW_1989_CONSTRAINED_ADVERSARIAL_AUDIT import optimal_sign_sequence_with_flip_budget
from QW_1993_GLOBAL_SCORE_COMPRESSION_HARD_LOCK_SEARCH import (
    TARGET_NULL_RANDOM_P90_MAX,
    TARGET_REAL_BLOCK_MIN,
    TARGET_REAL_IID_MIN,
    bootstrap_pass_rate,
    build_components,
)


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1997 = ROOT / "report_qw1997_fold2_guarded_seed_robust_hard_search.json"
OUT_JSON = ROOT / "report_qw1999_bounded_coupling_fold2_guarded_search.json"
OUT_MD = ROOT / "RAPORT_QW1999_BOUNDED_COUPLING_FOLD2_GUARDED_SEARCH.md"

# Fast stage budgets.
FAST_REAL_IID_BOOT = 60
FAST_NULL_TRIALS = 3
FAST_NULL_BOOT = 16
FAST_FOLD2_SEEDS = [290000, 290400, 290800]

# Full stage budgets.
FULL_REAL_IID_BOOT = 900
FULL_NULL_TRIALS = 14
FULL_NULL_BOOT = 70
FULL_SEEDS = [291000, 292000, 293000, 294000]
FULL_REAL_BLOCK_BOOT = 900
FULL_REAL_BLOCK_LEN = 10
FULL_ADV_FLIP_BUDGET = 4

SHORTLIST_SIZE = 12


def compress_score_power(s_raw: np.ndarray, gamma_c: float) -> np.ndarray:
    med = float(np.median(s_raw))
    std = max(float(np.std(s_raw)), 1e-12)
    z = (s_raw - med) / std
    zc = np.sign(z) * (np.abs(z) ** gamma_c)
    return med + std * zc


def bounded_coupling_term(t: np.ndarray, kappa_t: float) -> np.ndarray:
    t_std = max(float(np.std(t)), 1e-12)
    cap = float(kappa_t) * t_std
    return np.clip(t, -cap, cap)


def score_bounded(
    base: np.ndarray,
    c1: np.ndarray,
    c3: np.ndarray,
    c4: np.ndarray,
    xi1: float,
    xi3: float,
    xi4: float,
    gamma_c: float,
    kappa_t: float,
) -> np.ndarray:
    t = xi1 * c1 + xi3 * c3 + xi4 * c4
    t_eff = bounded_coupling_term(t, kappa_t=kappa_t)
    s_raw = base + t_eff
    return compress_score_power(s_raw, gamma_c=gamma_c)


def eval_fold_bounded(
    dff: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
    thr: Dict[str, float],
    xi1: float,
    xi3: float,
    xi4: float,
    gamma_c: float,
    kappa_t: float,
    real_iid_boot: int,
    null_trials: int,
    null_boot: int,
    seed: int,
) -> Dict[str, float]:
    base, c1, c3, c4, pairs = build_components(dff, kernel, params)
    s = score_bounded(
        base=base,
        c1=c1,
        c3=c3,
        c4=c4,
        xi1=xi1,
        xi3=xi3,
        xi4=xi4,
        gamma_c=gamma_c,
        kappa_t=kappa_t,
    )
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

        s_n = score_bounded(
            base=base,
            c1=c1n,
            c3=c3n,
            c4=c4n,
            xi1=xi1,
            xi3=xi3,
            xi4=xi4,
            gamma_c=gamma_c,
            kappa_t=kappa_t,
        )
        s_hl_n, s_hv_n, s_lv_n = s_n[pairs == 0], s_n[pairs == 1], s_n[pairs == 2]
        rb = bootstrap_pass_rate(s_hl_n, s_hv_n, s_lv_n, thr, n_boot=null_boot, seed=seed + 200 + t)
        rates.append(float(rb))

    arr = np.array(rates, dtype=float)
    return {
        "real_iid": float(real_iid),
        "null_random_mean": float(np.mean(arr)),
        "null_random_p90": float(np.quantile(arr, 0.90)),
    }


def constrained_adv_rate_bounded(
    dff: pd.DataFrame,
    kernel: Dict[str, float],
    params: Dict[str, float],
    thr: Dict[str, float],
    xi1: float,
    xi3: float,
    xi4: float,
    gamma_c: float,
    kappa_t: float,
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

    s_n = score_bounded(
        base=base,
        c1=c1n,
        c3=c3n,
        c4=c4n,
        xi1=xi1,
        xi3=xi3,
        xi4=xi4,
        gamma_c=gamma_c,
        kappa_t=kappa_t,
    )
    s_hl_n, s_hv_n, s_lv_n = s_n[pairs == 0], s_n[pairs == 1], s_n[pairs == 2]
    return float(bootstrap_pass_rate(s_hl_n, s_hv_n, s_lv_n, thr, n_boot=350, seed=seed))


def dedup_rows(rows: List[Dict[str, float]]) -> List[Dict[str, float]]:
    seen = set()
    out = []
    for r in rows:
        k = (
            round(float(r["xi1"]), 12),
            round(float(r["xi3"]), 12),
            round(float(r["xi4"]), 12),
            round(float(r["gamma_c"]), 6),
            round(float(r["kappa_t"]), 4),
        )
        if k in seen:
            continue
        seen.add(k)
        out.append({"xi1": k[0], "xi3": k[1], "xi4": k[2], "gamma_c": k[3], "kappa_t": k[4]})
    return out


def evaluate_candidate_fast(
    fold_dfs: List[pd.DataFrame],
    kernel: Dict[str, float],
    params: Dict[str, float],
    thr: Dict[str, float],
    xi1: float,
    xi3: float,
    xi4: float,
    gamma_c: float,
    kappa_t: float,
    seed_base: int,
) -> Dict[str, object]:
    fold_rows = []
    for f, dff in enumerate(fold_dfs):
        if f == 2:
            real_vals = []
            p90_vals = []
            for j, s in enumerate(FAST_FOLD2_SEEDS):
                row = eval_fold_bounded(
                    dff=dff,
                    kernel=kernel,
                    params=params,
                    thr=thr,
                    xi1=xi1,
                    xi3=xi3,
                    xi4=xi4,
                    gamma_c=gamma_c,
                    kappa_t=kappa_t,
                    real_iid_boot=FAST_REAL_IID_BOOT,
                    null_trials=FAST_NULL_TRIALS,
                    null_boot=FAST_NULL_BOOT,
                    seed=s + seed_base + j,
                )
                real_vals.append(float(row["real_iid"]))
                p90_vals.append(float(row["null_random_p90"]))
            fold_rows.append(
                {
                    "fold": f,
                    "real_iid_fast": float(np.mean(real_vals)),
                    "null_random_p90_fast": float(np.mean(p90_vals)),
                    "fold2_p90_std_fast": float(np.std(np.array(p90_vals, dtype=float))),
                }
            )
        else:
            row = eval_fold_bounded(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                gamma_c=gamma_c,
                kappa_t=kappa_t,
                real_iid_boot=FAST_REAL_IID_BOOT,
                null_trials=FAST_NULL_TRIALS,
                null_boot=FAST_NULL_BOOT,
                seed=seed_base + 10 * f,
            )
            fold_rows.append(
                {
                    "fold": f,
                    "real_iid_fast": float(row["real_iid"]),
                    "null_random_p90_fast": float(row["null_random_p90"]),
                    "fold2_p90_std_fast": 0.0,
                }
            )

    min_real = float(min(r["real_iid_fast"] for r in fold_rows))
    max_null = float(max(r["null_random_p90_fast"] for r in fold_rows))
    fold2_p90 = float([r["null_random_p90_fast"] for r in fold_rows if r["fold"] == 2][0])
    fold2_std = float([r["fold2_p90_std_fast"] for r in fold_rows if r["fold"] == 2][0])
    hard_margin = float(min(min_real - TARGET_REAL_IID_MIN, TARGET_NULL_RANDOM_P90_MAX - max_null))
    targeted = float(
        0.45 * (TARGET_NULL_RANDOM_P90_MAX - fold2_p90)
        + 0.25 * (TARGET_NULL_RANDOM_P90_MAX - max_null)
        + 0.20 * hard_margin
        - 0.10 * fold2_std
    )
    aux = float(np.mean([r["real_iid_fast"] for r in fold_rows]) - np.mean([r["null_random_p90_fast"] for r in fold_rows]))

    return {
        "fold_fast": fold_rows,
        "min_real_iid_fast": min_real,
        "max_null_random_p90_fast": max_null,
        "fold2_null_random_p90_fast": fold2_p90,
        "fold2_p90_std_fast": fold2_std,
        "hard_margin_fast": hard_margin,
        "targeted_score_fast": targeted,
        "aux_fast": aux,
    }


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    r1997 = json.loads(IN_QW1997.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]
    center = r1997["best"]

    rows = []
    d13 = np.array([-0.00004, -0.00002, 0.0, 0.00002, 0.00004], dtype=float)
    d4 = np.array([-0.00004, 0.0, 0.00004], dtype=float)
    dgc = np.array([-0.020, -0.010, 0.0, 0.010, 0.020], dtype=float)
    kappa_grid = np.array([1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.80], dtype=float)

    for a in d13:
        for b in d13:
            for d in d4:
                for g in dgc:
                    for kt in kappa_grid:
                        rows.append(
                            {
                                "xi1": float(np.clip(float(center["xi1"]) + a, 0.0001, 0.0032)),
                                "xi3": float(np.clip(float(center["xi3"]) + b, -0.0002, 0.0032)),
                                "xi4": float(np.clip(float(center["xi4"]) + d, -0.0032, 0.0032)),
                                "gamma_c": float(np.clip(float(center["gamma_c"]) + g, 0.75, 1.02)),
                                "kappa_t": float(kt),
                            }
                        )

    rng = np.random.default_rng(294000)
    for _ in range(180):
        rows.append(
            {
                "xi1": float(np.clip(rng.normal(float(center["xi1"]), 0.000022), 0.0001, 0.0032)),
                "xi3": float(np.clip(rng.normal(float(center["xi3"]), 0.000022), -0.0002, 0.0032)),
                "xi4": float(np.clip(rng.normal(float(center["xi4"]), 0.00003), -0.0032, 0.0032)),
                "gamma_c": float(np.clip(rng.normal(float(center["gamma_c"]), 0.014), 0.75, 1.02)),
                "kappa_t": float(np.clip(rng.normal(1.45, 0.22), 0.90, 2.20)),
            }
        )

    candidates = dedup_rows(rows)

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    fast_rows = []
    total = len(candidates)
    for i, c in enumerate(candidates):
        xi1 = float(c["xi1"])
        xi3 = float(c["xi3"])
        xi4 = float(c["xi4"])
        gc = float(c["gamma_c"])
        kt = float(c["kappa_t"])
        fr = evaluate_candidate_fast(
            fold_dfs=fold_dfs,
            kernel=kernel,
            params=params,
            thr=thr,
            xi1=xi1,
            xi3=xi3,
            xi4=xi4,
            gamma_c=gc,
            kappa_t=kt,
            seed_base=295000 + i * 20,
        )
        fast_rows.append({"xi1": xi1, "xi3": xi3, "xi4": xi4, "gamma_c": gc, "kappa_t": kt, **fr})
        if (i + 1) % 120 == 0:
            print(f"[QW-1999] fast search progress: {i + 1}/{total}", flush=True)

    fast_rows.sort(
        key=lambda x: (
            x["targeted_score_fast"],
            x["hard_margin_fast"],
            x["aux_fast"],
        ),
        reverse=True,
    )
    shortlist = fast_rows[:SHORTLIST_SIZE]
    print(f"[QW-1999] shortlist size: {len(shortlist)}", flush=True)

    final_rows = []
    for i, c in enumerate(shortlist):
        xi1 = float(c["xi1"])
        xi3 = float(c["xi3"])
        xi4 = float(c["xi4"])
        gc = float(c["gamma_c"])
        kt = float(c["kappa_t"])
        fold_results = []
        for f, dff in enumerate(fold_dfs):
            real_vals = []
            p90_vals = []
            mean_vals = []
            for j, s in enumerate(FULL_SEEDS):
                row = eval_fold_bounded(
                    dff=dff,
                    kernel=kernel,
                    params=params,
                    thr=thr,
                    xi1=xi1,
                    xi3=xi3,
                    xi4=xi4,
                    gamma_c=gc,
                    kappa_t=kt,
                    real_iid_boot=FULL_REAL_IID_BOOT,
                    null_trials=FULL_NULL_TRIALS,
                    null_boot=FULL_NULL_BOOT,
                    seed=s + 220 * f + j,
                )
                real_vals.append(float(row["real_iid"]))
                p90_vals.append(float(row["null_random_p90"]))
                mean_vals.append(float(row["null_random_mean"]))

            base_score, cc1, cc3, cc4, pairs = build_components(dff, kernel, params)
            s_real = score_bounded(
                base=base_score,
                c1=cc1,
                c3=cc3,
                c4=cc4,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                gamma_c=gc,
                kappa_t=kt,
            )
            s_hl, s_hv, s_lv = s_real[pairs == 0], s_real[pairs == 1], s_real[pairs == 2]
            real_block = bootstrap_pass_rate_block(
                s_hl=s_hl,
                s_hv=s_hv,
                s_lv=s_lv,
                thr=thr,
                n_boot=FULL_REAL_BLOCK_BOOT,
                seed=299000 + i * 50 + f,
                block_len=FULL_REAL_BLOCK_LEN,
            )
            adv = constrained_adv_rate_bounded(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                gamma_c=gc,
                kappa_t=kt,
                max_flips=FULL_ADV_FLIP_BUDGET,
                seed=300000 + i * 50 + f,
            )

            arr_r = np.array(real_vals, dtype=float)
            arr_p90 = np.array(p90_vals, dtype=float)
            arr_m = np.array(mean_vals, dtype=float)
            fold_results.append(
                {
                    "fold": f,
                    "real_iid_mean_full": float(np.mean(arr_r)),
                    "real_iid_p10_full": float(np.quantile(arr_r, 0.10)),
                    "real_iid_min_full": float(np.min(arr_r)),
                    "null_random_mean_mean_full": float(np.mean(arr_m)),
                    "null_random_p90_mean_full": float(np.mean(arr_p90)),
                    "null_random_p90_p90_full": float(np.quantile(arr_p90, 0.90)),
                    "null_random_p90_max_full": float(np.max(arr_p90)),
                    "real_block_full": float(real_block),
                    "adv_constrained_full": float(adv),
                }
            )

        min_real_iid_mean = float(min(r["real_iid_mean_full"] for r in fold_results))
        min_real_iid_p10 = float(min(r["real_iid_p10_full"] for r in fold_results))
        min_real_block = float(min(r["real_block_full"] for r in fold_results))
        max_null_p90_mean = float(max(r["null_random_p90_mean_full"] for r in fold_results))
        max_null_p90_p90 = float(max(r["null_random_p90_p90_full"] for r in fold_results))
        max_adv = float(max(r["adv_constrained_full"] for r in fold_results))
        fold2_mean = float([r["null_random_p90_mean_full"] for r in fold_results if r["fold"] == 2][0])
        fold2_p90 = float([r["null_random_p90_p90_full"] for r in fold_results if r["fold"] == 2][0])

        hard_margin = float(min(min_real_iid_mean - TARGET_REAL_IID_MIN, TARGET_NULL_RANDOM_P90_MAX - max_null_p90_mean))
        hard_pass = bool(
            min_real_iid_mean >= TARGET_REAL_IID_MIN
            and min_real_iid_p10 >= 0.94
            and min_real_block >= TARGET_REAL_BLOCK_MIN
            and max_null_p90_mean <= TARGET_NULL_RANDOM_P90_MAX
            and max_null_p90_p90 <= 0.45
            and fold2_p90 <= 0.45
            and max_adv <= 0.45
        )
        targeted = float(
            0.45 * (TARGET_NULL_RANDOM_P90_MAX - fold2_mean)
            + 0.30 * (TARGET_NULL_RANDOM_P90_MAX - max_null_p90_mean)
            + 0.15 * (0.45 - max_null_p90_p90)
            + 0.10 * hard_margin
        )
        aux = float(
            np.mean([r["real_iid_mean_full"] for r in fold_results])
            - np.mean([r["null_random_p90_mean_full"] for r in fold_results])
        )

        final_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "xi4": xi4,
                "gamma_c": gc,
                "kappa_t": kt,
                "hard_pass": hard_pass,
                "min_real_iid_mean_full": min_real_iid_mean,
                "min_real_iid_p10_full": min_real_iid_p10,
                "min_real_block_full": min_real_block,
                "max_null_random_p90_mean_full": max_null_p90_mean,
                "max_null_random_p90_p90_full": max_null_p90_p90,
                "max_adv_constrained_full": max_adv,
                "fold2_null_random_p90_mean_full": fold2_mean,
                "fold2_null_random_p90_p90_full": fold2_p90,
                "hard_margin_full": hard_margin,
                "targeted_score_full": targeted,
                "aux_full": aux,
                "fold_results": fold_results,
            }
        )
        print(f"[QW-1999] full eval done: {i + 1}/{len(shortlist)}", flush=True)

    final_rows.sort(
        key=lambda x: (
            int(x["hard_pass"]),
            x["targeted_score_full"],
            x["hard_margin_full"],
            x["aux_full"],
        ),
        reverse=True,
    )
    best = final_rows[0]
    verdict = "BOUNDED_COUPLING_ROBUST_HARD_PASS" if best["hard_pass"] else "BOUNDED_COUPLING_ROBUST_HARD_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1997.name],
        "search_config": {
            "total_candidates": total,
            "shortlist_size": SHORTLIST_SIZE,
            "target_real_iid_min": TARGET_REAL_IID_MIN,
            "target_null_random_p90_max": TARGET_NULL_RANDOM_P90_MAX,
            "target_real_block_min": TARGET_REAL_BLOCK_MIN,
            "fast_fold2_seeds": FAST_FOLD2_SEEDS,
            "full_seeds": FULL_SEEDS,
            "full_adv_flip_budget": FULL_ADV_FLIP_BUDGET,
        },
        "best": best,
        "top12": final_rows,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1999: BOUNDED COUPLING FOLD-2 GUARDED SEARCH",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Best Candidate",
        f"- xi1/xi3/xi4: {best['xi1']:.6f}/{best['xi3']:.6f}/{best['xi4']:.6f}",
        f"- gamma_c/kappa_t: {best['gamma_c']:.4f}/{best['kappa_t']:.3f}",
        f"- hard_pass: {best['hard_pass']}",
        f"- min_real_iid_mean_full: {100.0 * best['min_real_iid_mean_full']:.2f}%",
        f"- min_real_iid_p10_full: {100.0 * best['min_real_iid_p10_full']:.2f}%",
        f"- min_real_block_full: {100.0 * best['min_real_block_full']:.2f}%",
        f"- max_null_random_p90_mean_full: {100.0 * best['max_null_random_p90_mean_full']:.2f}%",
        f"- max_null_random_p90_p90_full: {100.0 * best['max_null_random_p90_p90_full']:.2f}%",
        f"- fold2 p90 mean/p90: {100.0 * best['fold2_null_random_p90_mean_full']:.2f}% / {100.0 * best['fold2_null_random_p90_p90_full']:.2f}%",
        f"- max_adv_constrained_full: {100.0 * best['max_adv_constrained_full']:.2f}%",
        f"- hard_margin_full: {100.0 * best['hard_margin_full']:.2f} pp",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1999] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1999] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1999] verdict={verdict}")


if __name__ == "__main__":
    main()
