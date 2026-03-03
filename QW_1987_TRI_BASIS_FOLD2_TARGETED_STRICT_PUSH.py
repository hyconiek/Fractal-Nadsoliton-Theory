#!/usr/bin/env python3
"""
QW-1987: Tri-basis strict push targeted at fold-2 null leakage.
No fold-specific retune; single global operator (xi1, xi3, xi4).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import ROOT
from QW_1986_TRI_BASIS_STRICT_5OF5_ATTEMPT import eval_on_fold


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1986 = ROOT / "report_qw1986_tri_basis_strict_5of5_attempt.json"
OUT_JSON = ROOT / "report_qw1987_tri_basis_fold2_targeted_strict_push.json"
OUT_MD = ROOT / "RAPORT_QW1987_TRI_BASIS_FOLD2_TARGETED_STRICT_PUSH.md"

REAL_PASS_MIN = 0.85
NULL_P90_PASS_MAX = 0.40

FAST_REAL_BOOT = 100
FAST_NULL_TRIALS = 6
FAST_NULL_BOOT = 30

FULL_REAL_BOOT = 1200
FULL_NULL_TRIALS = 16
FULL_NULL_BOOT = 80
SHORTLIST_SIZE = 16


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
    r1986 = json.loads(IN_QW1986.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    xi1_c = float(r1986["best"]["xi1"])
    xi3_c = float(r1986["best"]["xi3"])
    xi4_c = float(r1986["best"]["xi4"])

    d13 = np.array([-0.00006, -0.00003, -0.00001, 0.0, 0.00001, 0.00003, 0.00006], dtype=float)
    d4 = np.array([-0.00025, -0.00012, -0.00006, 0.0, 0.00006, 0.00012, 0.00025], dtype=float)

    rows = []
    for a in d13:
        for b in d13:
            for c in d4:
                rows.append(
                    {
                        "xi1": float(np.clip(xi1_c + a, 0.0001, 0.0030)),
                        "xi3": float(np.clip(xi3_c + b, -0.0002, 0.0030)),
                        "xi4": float(np.clip(xi4_c + c, -0.0030, 0.0030)),
                    }
                )
    rng = np.random.default_rng(170000)
    for _ in range(220):
        rows.append(
            {
                "xi1": float(np.clip(rng.normal(xi1_c, 0.00005), 0.0001, 0.0030)),
                "xi3": float(np.clip(rng.normal(xi3_c, 0.00005), -0.0002, 0.0030)),
                "xi4": float(np.clip(rng.normal(xi4_c, 0.00015), -0.0030, 0.0030)),
            }
        )
    candidates = unique_rows(rows)

    fast_rows = []
    total = len(candidates)
    for i, c in enumerate(candidates):
        xi1 = float(c["xi1"])
        xi3 = float(c["xi3"])
        xi4 = float(c["xi4"])
        folds = []
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
                seed=171000 + i * 10 + f,
            )
            folds.append(
                {
                    "fold": f,
                    "real_fast": row["real_rate"],
                    "null_p90_fast": row["null_p90"],
                    "null_mean_fast": row["null_mean"],
                    "det_flags": row["det_flags"],
                }
            )

        min_real = float(min(r["real_fast"] for r in folds))
        max_null = float(max(r["null_p90_fast"] for r in folds))
        fold2_null = float([r["null_p90_fast"] for r in folds if r["fold"] == 2][0])
        strict_margin = float(min(min_real - REAL_PASS_MIN, NULL_P90_PASS_MAX - max_null))
        strict_fast = bool(strict_margin >= 0.0)
        # Targeted score: worst-fold strict margin + extra pressure on fold 2.
        targeted_score = float(strict_margin + 0.5 * (NULL_P90_PASS_MAX - fold2_null))
        aux = float(np.mean([r["real_fast"] for r in folds]) - np.mean([r["null_p90_fast"] for r in folds]))
        fast_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "xi4": xi4,
                "fold_fast": folds,
                "min_real_fast": min_real,
                "max_null_p90_fast": max_null,
                "fold2_null_p90_fast": fold2_null,
                "strict_margin_fast": strict_margin,
                "strict_fast": strict_fast,
                "targeted_score_fast": targeted_score,
                "aux_score_fast": aux,
            }
        )
        if (i + 1) % 80 == 0:
            print(f"[QW-1987] fast search progress: {i + 1}/{total}", flush=True)

    fast_rows.sort(
        key=lambda x: (
            int(x["strict_fast"]),
            x["targeted_score_fast"],
            x["aux_score_fast"],
        ),
        reverse=True,
    )
    shortlist = fast_rows[:SHORTLIST_SIZE]
    print(f"[QW-1987] shortlist size: {len(shortlist)}", flush=True)

    final_rows = []
    for i, c in enumerate(shortlist):
        xi1 = float(c["xi1"])
        xi3 = float(c["xi3"])
        xi4 = float(c["xi4"])
        folds = []
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
                seed=175000 + i * 10 + f,
            )
            fold_pass = bool(row["real_rate"] >= REAL_PASS_MIN and row["null_p90"] <= NULL_P90_PASS_MAX)
            folds.append(
                {
                    "fold": f,
                    "real_full": row["real_rate"],
                    "null_mean_full": row["null_mean"],
                    "null_p90_full": row["null_p90"],
                    "det_flags": row["det_flags"],
                    "fold_pass": fold_pass,
                }
            )
        pass_count = int(sum(int(r["fold_pass"]) for r in folds))
        min_real = float(min(r["real_full"] for r in folds))
        max_null = float(max(r["null_p90_full"] for r in folds))
        fold2_null = float([r["null_p90_full"] for r in folds if r["fold"] == 2][0])
        strict_margin = float(min(min_real - REAL_PASS_MIN, NULL_P90_PASS_MAX - max_null))
        targeted_score = float(strict_margin + 0.5 * (NULL_P90_PASS_MAX - fold2_null))
        aux = float(np.mean([r["real_full"] for r in folds]) - np.mean([r["null_p90_full"] for r in folds]))
        final_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "xi4": xi4,
                "pass_count": pass_count,
                "min_real_full": min_real,
                "max_null_p90_full": max_null,
                "fold2_null_p90_full": fold2_null,
                "strict_margin_full": strict_margin,
                "targeted_score_full": targeted_score,
                "aux_score_full": aux,
                "fold_results": folds,
            }
        )

    final_rows.sort(
        key=lambda x: (
            x["pass_count"],
            x["targeted_score_full"],
            x["aux_score_full"],
        ),
        reverse=True,
    )
    best = final_rows[0]

    verdict = (
        "TRI_BASIS_TARGETED_STRICT_5OF5_PASS"
        if best["pass_count"] == 5
        else "TRI_BASIS_TARGETED_NEAR_MISS"
        if best["pass_count"] == 4
        else "TRI_BASIS_TARGETED_FAIL"
    )
    required = (
        "LOCK_OPERATOR_AND_RUN_HARDER_INDEPENDENT_CONFIRMATORY_SUITE"
        if best["pass_count"] == 5
        else "BORDERLINE_REMAINS_NEED_STRONGER_INVARIANT_OR_REFORMULATE_TEST_STATISTIC"
        if best["pass_count"] == 4
        else "CURRENT_TRI_BASIS_UNSUITABLE_FOR_STRICT_5OF5"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1986.name],
        "search_config": {
            "grid_points": len(d13) * len(d13) * len(d4),
            "random_points": 220,
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
        "top16": final_rows,
        "verdict": verdict,
        "required_next_step": required,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1987: TRI-BASIS FOLD-2 TARGETED STRICT PUSH",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Best Candidate",
        f"- xi1/xi3/xi4: {best['xi1']:.6f}/{best['xi3']:.6f}/{best['xi4']:.6f}",
        f"- pass_count: {best['pass_count']}/5",
        f"- min_real_full: {100.0 * best['min_real_full']:.2f}%",
        f"- max_null_p90_full: {100.0 * best['max_null_p90_full']:.2f}%",
        f"- fold2_null_p90_full: {100.0 * best['fold2_null_p90_full']:.2f}%",
        f"- strict_margin_full: {100.0 * best['strict_margin_full']:.2f} pp",
        "",
        "## Required Next Step",
        f"- {required}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1987] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1987] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1987] verdict={verdict} pass_count={best['pass_count']}/5")


if __name__ == "__main__":
    main()

