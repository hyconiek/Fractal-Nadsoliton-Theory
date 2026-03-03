#!/usr/bin/env python3
"""
QW-1998: Deep audit for QW-1997 fold-2 guarded robust hard-pass candidate.
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
    TARGET_REAL_IID_MIN,
    constrained_adv_rate_compressed,
    eval_fold,
)


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1997 = ROOT / "report_qw1997_fold2_guarded_seed_robust_hard_search.json"
OUT_JSON = ROOT / "report_qw1998_fold2_guarded_deep_audit.json"
OUT_MD = ROOT / "RAPORT_QW1998_FOLD2_GUARDED_DEEP_AUDIT.md"

SEEDS = [280000, 281000, 282000, 283000, 284000, 285000]
REAL_IID_BOOT = 1200
NULL_TRIALS = 40
NULL_BOOT = 120
BLOCK_BOOT = 1200
BLOCK_LEN = 10
ADV_FLIP_BUDGET = 4


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    r1997 = json.loads(IN_QW1997.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]
    b = r1997["best"]
    xi1 = float(b["xi1"])
    xi3 = float(b["xi3"])
    xi4 = float(b["xi4"])
    gamma_c = float(b["gamma_c"])

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    # Needed for block and constrained-adversarial.
    from QW_1993_GLOBAL_SCORE_COMPRESSION_HARD_LOCK_SEARCH import build_components, compress_score

    fold_rows = []
    for f, dff in enumerate(fold_dfs):
        real_vals = []
        null_mean_vals = []
        null_p90_vals = []
        for j, s in enumerate(SEEDS):
            row = eval_fold(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                gamma_c=gamma_c,
                real_iid_boot=REAL_IID_BOOT,
                null_trials=NULL_TRIALS,
                null_boot=NULL_BOOT,
                seed=s + 100 * f + j,
            )
            real_vals.append(float(row["real_iid"]))
            null_mean_vals.append(float(row["null_random_mean"]))
            null_p90_vals.append(float(row["null_random_p90"]))

        base_score, c1, c3, c4, pairs = build_components(dff, kernel, params)
        s_raw = base_score + xi1 * c1 + xi3 * c3 + xi4 * c4
        s = compress_score(s_raw, gamma_c)
        s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]
        real_block = bootstrap_pass_rate_block(
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            thr=thr,
            n_boot=BLOCK_BOOT,
            seed=286000 + f,
            block_len=BLOCK_LEN,
        )
        adv = constrained_adv_rate_compressed(
            dff=dff,
            kernel=kernel,
            params=params,
            thr=thr,
            xi1=xi1,
            xi3=xi3,
            xi4=xi4,
            gamma_c=gamma_c,
            max_flips=ADV_FLIP_BUDGET,
            seed=287000 + f,
        )

        arr_r = np.array(real_vals, dtype=float)
        arr_nm = np.array(null_mean_vals, dtype=float)
        arr_np = np.array(null_p90_vals, dtype=float)
        fold_rows.append(
            {
                "fold": f,
                "real_iid_mean": float(np.mean(arr_r)),
                "real_iid_p10": float(np.quantile(arr_r, 0.10)),
                "real_iid_min": float(np.min(arr_r)),
                "null_random_mean_mean": float(np.mean(arr_nm)),
                "null_random_p90_mean": float(np.mean(arr_np)),
                "null_random_p90_p90": float(np.quantile(arr_np, 0.90)),
                "null_random_p90_max": float(np.max(arr_np)),
                "real_block": float(real_block),
                "adv_constrained": float(adv),
            }
        )
        print(f"[QW-1998] fold {f} done", flush=True)

    agg = {
        "real_iid_mean_min": float(min(r["real_iid_mean"] for r in fold_rows)),
        "real_iid_p10_min": float(min(r["real_iid_p10"] for r in fold_rows)),
        "real_iid_min_min": float(min(r["real_iid_min"] for r in fold_rows)),
        "null_random_p90_mean_max": float(max(r["null_random_p90_mean"] for r in fold_rows)),
        "null_random_p90_p90_max": float(max(r["null_random_p90_p90"] for r in fold_rows)),
        "null_random_p90_max_max": float(max(r["null_random_p90_max"] for r in fold_rows)),
        "real_block_min": float(min(r["real_block"] for r in fold_rows)),
        "adv_constrained_max": float(max(r["adv_constrained"] for r in fold_rows)),
    }
    hard_pass = bool(
        agg["real_iid_mean_min"] >= TARGET_REAL_IID_MIN
        and agg["real_iid_p10_min"] >= 0.94
        and agg["null_random_p90_mean_max"] <= TARGET_NULL_RANDOM_P90_MAX
        and agg["null_random_p90_p90_max"] <= 0.45
        and agg["real_block_min"] >= 0.90
        and agg["adv_constrained_max"] <= 0.45
    )
    near_pass = bool(
        agg["real_iid_mean_min"] >= 0.948
        and agg["real_iid_p10_min"] >= 0.938
        and agg["null_random_p90_mean_max"] <= 0.42
        and agg["null_random_p90_p90_max"] <= 0.47
        and agg["real_block_min"] >= 0.90
        and agg["adv_constrained_max"] <= 0.45
    )
    verdict = (
        "FOLD2_GUARDED_DEEP_AUDIT_PASS"
        if hard_pass
        else "FOLD2_GUARDED_DEEP_AUDIT_NEAR_PASS"
        if near_pass
        else "FOLD2_GUARDED_DEEP_AUDIT_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1997.name],
        "candidate": {"xi1": xi1, "xi3": xi3, "xi4": xi4, "gamma_c": gamma_c},
        "config": {
            "seeds": SEEDS,
            "real_iid_boot": REAL_IID_BOOT,
            "null_trials": NULL_TRIALS,
            "null_boot": NULL_BOOT,
            "block_boot": BLOCK_BOOT,
            "adv_flip_budget": ADV_FLIP_BUDGET,
        },
        "fold_results": fold_rows,
        "aggregate": agg,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1998: FOLD-2 GUARDED DEEP AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Aggregate",
        f"- real_iid_mean_min: {100.0 * agg['real_iid_mean_min']:.2f}%",
        f"- real_iid_p10_min: {100.0 * agg['real_iid_p10_min']:.2f}%",
        f"- real_iid_min_min: {100.0 * agg['real_iid_min_min']:.2f}%",
        f"- null_random_p90_mean_max: {100.0 * agg['null_random_p90_mean_max']:.2f}%",
        f"- null_random_p90_p90_max: {100.0 * agg['null_random_p90_p90_max']:.2f}%",
        f"- null_random_p90_max_max: {100.0 * agg['null_random_p90_max_max']:.2f}%",
        f"- real_block_min: {100.0 * agg['real_block_min']:.2f}%",
        f"- adv_constrained_max: {100.0 * agg['adv_constrained_max']:.2f}%",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1998] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1998] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1998] verdict={verdict}")


if __name__ == "__main__":
    main()

