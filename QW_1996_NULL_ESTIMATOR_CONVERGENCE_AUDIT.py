#!/usr/bin/env python3
"""
QW-1996: Null-estimator convergence audit for QW-1994 candidate.
Checks whether hard-pass fragility is estimator noise or true instability.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import ROOT
from QW_1993_GLOBAL_SCORE_COMPRESSION_HARD_LOCK_SEARCH import eval_fold


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1994 = ROOT / "report_qw1994_compression_micro_local_hard_push.json"
OUT_JSON = ROOT / "report_qw1996_null_estimator_convergence_audit.json"
OUT_MD = ROOT / "RAPORT_QW1996_NULL_ESTIMATOR_CONVERGENCE_AUDIT.md"

SEEDS = [260000, 261000, 262000]
REAL_IID_BOOT = 1000
NULL_TRIALS = 48
NULL_BOOT = 100


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    r1994 = json.loads(IN_QW1994.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]
    b = r1994["best"]
    xi1 = float(b["xi1"])
    xi3 = float(b["xi3"])
    xi4 = float(b["xi4"])
    gamma_c = float(b["gamma_c"])

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    fold_rows = []
    for f, dff in enumerate(fold_dfs):
        real_vals = []
        p90_vals = []
        mean_vals = []
        for i, s in enumerate(SEEDS):
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
                seed=s + 100 * f + i,
            )
            real_vals.append(float(row["real_iid"]))
            p90_vals.append(float(row["null_random_p90"]))
            mean_vals.append(float(row["null_random_mean"]))

        arr_r = np.array(real_vals, dtype=float)
        arr_p90 = np.array(p90_vals, dtype=float)
        arr_m = np.array(mean_vals, dtype=float)
        fold_rows.append(
            {
                "fold": f,
                "real_iid_mean": float(np.mean(arr_r)),
                "real_iid_min": float(np.min(arr_r)),
                "null_random_mean_mean": float(np.mean(arr_m)),
                "null_random_p90_mean": float(np.mean(arr_p90)),
                "null_random_p90_min": float(np.min(arr_p90)),
                "null_random_p90_max": float(np.max(arr_p90)),
                "null_random_p90_std": float(np.std(arr_p90)),
            }
        )
        print(f"[QW-1996] fold {f} done", flush=True)

    agg = {
        "real_iid_min_overall": float(min(r["real_iid_min"] for r in fold_rows)),
        "null_random_p90_mean_max": float(max(r["null_random_p90_mean"] for r in fold_rows)),
        "null_random_p90_max_max": float(max(r["null_random_p90_max"] for r in fold_rows)),
        "null_random_p90_std_max": float(max(r["null_random_p90_std"] for r in fold_rows)),
    }
    verdict = (
        "NULL_ESTIMATOR_CONVERGED_PASSLIKE"
        if agg["null_random_p90_mean_max"] <= 0.40 and agg["null_random_p90_max_max"] <= 0.45
        else "NULL_ESTIMATOR_CONVERGED_FAILLIKE"
        if agg["null_random_p90_std_max"] <= 0.05
        else "NULL_ESTIMATOR_STILL_HIGH_VARIANCE"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1994.name],
        "candidate": {"xi1": xi1, "xi3": xi3, "xi4": xi4, "gamma_c": gamma_c},
        "config": {
            "seeds": SEEDS,
            "real_iid_boot": REAL_IID_BOOT,
            "null_trials": NULL_TRIALS,
            "null_boot": NULL_BOOT,
        },
        "fold_results": fold_rows,
        "aggregate": agg,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1996: NULL ESTIMATOR CONVERGENCE AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Aggregate",
        f"- real_iid_min_overall: {100.0 * agg['real_iid_min_overall']:.2f}%",
        f"- null_random_p90_mean_max: {100.0 * agg['null_random_p90_mean_max']:.2f}%",
        f"- null_random_p90_max_max: {100.0 * agg['null_random_p90_max_max']:.2f}%",
        f"- null_random_p90_std_max: {100.0 * agg['null_random_p90_std_max']:.2f}%",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1996] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1996] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1996] verdict={verdict}")


if __name__ == "__main__":
    main()

