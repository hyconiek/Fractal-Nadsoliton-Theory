#!/usr/bin/env python3
"""
QW-1997: Fold-2 guarded seed-robust hard-lock search.
Searches a single frozen operator with explicit robustness against fold-2 null instability.
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
from QW_1993_GLOBAL_SCORE_COMPRESSION_HARD_LOCK_SEARCH import (
    TARGET_NULL_RANDOM_P90_MAX,
    TARGET_REAL_BLOCK_MIN,
    TARGET_REAL_IID_MIN,
    constrained_adv_rate_compressed,
    eval_fold,
)


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1994 = ROOT / "report_qw1994_compression_micro_local_hard_push.json"
OUT_JSON = ROOT / "report_qw1997_fold2_guarded_seed_robust_hard_search.json"
OUT_MD = ROOT / "RAPORT_QW1997_FOLD2_GUARDED_SEED_ROBUST_HARD_SEARCH.md"

# Fast stage budgets (cheap, but with fold-2 seed averaging).
FAST_REAL_IID_BOOT = 80
FAST_NULL_TRIALS = 4
FAST_NULL_BOOT = 20
FAST_FOLD2_SEEDS = [270000, 270500]

# Full stage budgets (seed-robust).
FULL_REAL_IID_BOOT = 1000
FULL_NULL_TRIALS = 16
FULL_NULL_BOOT = 80
FULL_SEEDS = [271000, 272000, 273000]
FULL_REAL_BLOCK_BOOT = 850
FULL_REAL_BLOCK_LEN = 10
FULL_ADV_FLIP_BUDGET = 4

SHORTLIST_SIZE = 10


def dedup_rows(rows: List[Dict[str, float]]) -> List[Dict[str, float]]:
    seen = set()
    out = []
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
        out.append({"xi1": k[0], "xi3": k[1], "xi4": k[2], "gamma_c": k[3]})
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
    seed_base: int,
) -> Dict[str, object]:
    fold_rows = []
    for f, dff in enumerate(fold_dfs):
        if f == 2:
            # Explicit guard: fold-2 evaluated on multiple seeds even in fast stage.
            real_vals = []
            p90_vals = []
            for j, s in enumerate(FAST_FOLD2_SEEDS):
                row = eval_fold(
                    dff=dff,
                    kernel=kernel,
                    params=params,
                    thr=thr,
                    xi1=xi1,
                    xi3=xi3,
                    xi4=xi4,
                    gamma_c=gamma_c,
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
            row = eval_fold(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                gamma_c=gamma_c,
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
    # Explicit fold-2 guard in objective: keep fold-2 p90 low and stable.
    targeted = float(
        0.60 * (TARGET_NULL_RANDOM_P90_MAX - fold2_p90)
        + 0.30 * hard_margin
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
    r1994 = json.loads(IN_QW1994.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]
    center = r1994["best"]

    rows = []
    d13 = np.array([-0.00005, -0.00003, -0.00001, 0.0, 0.00001, 0.00003, 0.00005], dtype=float)
    d4 = np.array([-0.00006, -0.00003, 0.0, 0.00003, 0.00006], dtype=float)
    dgc = np.array([-0.020, -0.010, -0.005, 0.0, 0.005, 0.010, 0.020], dtype=float)
    for a in d13:
        for b in d13:
            for d in d4:
                for g in dgc:
                    rows.append(
                        {
                            "xi1": float(np.clip(float(center["xi1"]) + a, 0.0001, 0.0032)),
                            "xi3": float(np.clip(float(center["xi3"]) + b, -0.0002, 0.0032)),
                            "xi4": float(np.clip(float(center["xi4"]) + d, -0.0032, 0.0032)),
                            "gamma_c": float(np.clip(float(center["gamma_c"]) + g, 0.78, 1.02)),
                        }
                    )

    rng = np.random.default_rng(275000)
    for _ in range(180):
        rows.append(
            {
                "xi1": float(np.clip(rng.normal(float(center["xi1"]), 0.00002), 0.0001, 0.0032)),
                "xi3": float(np.clip(rng.normal(float(center["xi3"]), 0.00002), -0.0002, 0.0032)),
                "xi4": float(np.clip(rng.normal(float(center["xi4"]), 0.00003), -0.0032, 0.0032)),
                "gamma_c": float(np.clip(rng.normal(float(center["gamma_c"]), 0.012), 0.78, 1.02)),
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
        xi1, xi3, xi4, gc = float(c["xi1"]), float(c["xi3"]), float(c["xi4"]), float(c["gamma_c"])
        fr = evaluate_candidate_fast(
            fold_dfs=fold_dfs,
            kernel=kernel,
            params=params,
            thr=thr,
            xi1=xi1,
            xi3=xi3,
            xi4=xi4,
            gamma_c=gc,
            seed_base=276000 + i * 20,
        )
        fast_rows.append({"xi1": xi1, "xi3": xi3, "xi4": xi4, "gamma_c": gc, **fr})
        if (i + 1) % 100 == 0:
            print(f"[QW-1997] fast search progress: {i + 1}/{total}", flush=True)

    fast_rows.sort(
        key=lambda x: (
            x["targeted_score_fast"],
            x["hard_margin_fast"],
            x["aux_fast"],
        ),
        reverse=True,
    )
    shortlist = fast_rows[:SHORTLIST_SIZE]
    print(f"[QW-1997] shortlist size: {len(shortlist)}", flush=True)

    final_rows = []
    for i, c in enumerate(shortlist):
        xi1, xi3, xi4, gc = float(c["xi1"]), float(c["xi3"]), float(c["xi4"]), float(c["gamma_c"])
        fold_results = []
        for f, dff in enumerate(fold_dfs):
            real_vals = []
            p90_vals = []
            mean_vals = []
            for j, s in enumerate(FULL_SEEDS):
                row = eval_fold(
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
                    seed=s + 200 * f + j,
                )
                real_vals.append(float(row["real_iid"]))
                p90_vals.append(float(row["null_random_p90"]))
                mean_vals.append(float(row["null_random_mean"]))

            # block + constrained adversarial on this fold (single high-budget eval)
            from QW_1993_GLOBAL_SCORE_COMPRESSION_HARD_LOCK_SEARCH import build_components, compress_score

            base_score, cc1, cc3, cc4, pairs = build_components(dff, kernel, params)
            s_raw = base_score + xi1 * cc1 + xi3 * cc3 + xi4 * cc4
            s = compress_score(s_raw, gc)
            s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]
            real_block = bootstrap_pass_rate_block(
                s_hl=s_hl,
                s_hv=s_hv,
                s_lv=s_lv,
                thr=thr,
                n_boot=FULL_REAL_BLOCK_BOOT,
                seed=279000 + i * 50 + f,
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
                seed=280000 + i * 50 + f,
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
            and min_real_block >= TARGET_REAL_BLOCK_MIN
            and max_null_p90_mean <= TARGET_NULL_RANDOM_P90_MAX
            and max_adv <= 0.45
            and fold2_p90 <= 0.45
            and min_real_iid_p10 >= 0.94
        )
        targeted = float(
            0.55 * (TARGET_NULL_RANDOM_P90_MAX - fold2_mean)
            + 0.25 * hard_margin
            + 0.20 * (TARGET_NULL_RANDOM_P90_MAX - max_null_p90_mean)
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
    verdict = "FOLD2_GUARDED_ROBUST_HARD_PASS" if best["hard_pass"] else "FOLD2_GUARDED_ROBUST_HARD_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1994.name],
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
        "top10": final_rows,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1997: FOLD-2 GUARDED SEED-ROBUST HARD SEARCH",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Best Candidate",
        f"- xi1/xi3/xi4/gamma_c: {best['xi1']:.6f}/{best['xi3']:.6f}/{best['xi4']:.6f}/{best['gamma_c']:.3f}",
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

    print(f"[QW-1997] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1997] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1997] verdict={verdict}")


if __name__ == "__main__":
    main()

