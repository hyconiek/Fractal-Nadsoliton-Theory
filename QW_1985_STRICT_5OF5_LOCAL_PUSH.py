#!/usr/bin/env python3
"""
QW-1985: Strict 5/5 local push around the best QW-1984 candidate.
Single frozen operator only (no fold retune).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import ROOT
from QW_1984_FOLD_MINMAX_NULLP90_REFINEMENT import eval_candidate_on_fold


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1984 = ROOT / "report_qw1984_fold_minmax_nullp90_refinement.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
OUT_JSON = ROOT / "report_qw1985_strict_5of5_local_push.json"
OUT_MD = ROOT / "RAPORT_QW1985_STRICT_5OF5_LOCAL_PUSH.md"

REAL_PASS_MIN = 0.85
NULL_P90_PASS_MAX = 0.40

FAST_REAL_BOOT = 100
FAST_NULL_TRIALS = 6
FAST_NULL_BOOT = 30

FULL_REAL_BOOT = 1000
FULL_NULL_TRIALS = 14
FULL_NULL_BOOT = 70

SHORTLIST_SIZE = 12


def unique_candidates(rows: List[Dict[str, float]]) -> List[Dict[str, float]]:
    seen = set()
    out = []
    for r in rows:
        k = (round(float(r["xi1"]), 12), round(float(r["xi3"]), 12))
        if k in seen:
            continue
        seen.add(k)
        out.append({"xi1": float(k[0]), "xi3": float(k[1])})
    return out


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1984 = json.loads(IN_QW1984.read_text(encoding="utf-8"))
    thr = json.loads(IN_QW1969.read_text(encoding="utf-8"))["thresholds"]
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    xi1_c = float(r1984["best"]["xi1"])
    xi3_c = float(r1984["best"]["xi3"])

    deltas = np.array(
        [
            -0.00024,
            -0.00018,
            -0.00012,
            -0.00008,
            -0.00004,
            -0.00002,
            0.0,
            0.00002,
            0.00004,
            0.00008,
            0.00012,
            0.00018,
            0.00024,
        ],
        dtype=float,
    )
    rows = []
    for d1 in deltas:
        for d3 in deltas:
            rows.append(
                {
                    "xi1": float(np.clip(xi1_c + d1, 0.0001, 0.0026)),
                    "xi3": float(np.clip(xi3_c + d3, -0.0002, 0.0025)),
                }
            )
    rng = np.random.default_rng(150000)
    for _ in range(120):
        rows.append(
            {
                "xi1": float(np.clip(rng.normal(xi1_c, 0.00010), 0.0001, 0.0026)),
                "xi3": float(np.clip(rng.normal(xi3_c, 0.00010), -0.0002, 0.0025)),
            }
        )
    candidates = unique_candidates(rows)

    fast_rows = []
    total = len(candidates)
    for i, c in enumerate(candidates):
        xi1 = c["xi1"]
        xi3 = c["xi3"]
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
                seed=151000 + i * 10 + f,
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
        max_null = float(max(r["null_p90_fast"] for r in fold_rows))
        strict_margin = float(min(min_real - REAL_PASS_MIN, NULL_P90_PASS_MAX - max_null))
        strict_fast = bool(strict_margin >= 0.0)
        aux = float(np.mean([r["real_fast"] for r in fold_rows]) - np.mean([r["null_p90_fast"] for r in fold_rows]))
        fast_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "fold_fast": fold_rows,
                "min_real_fast": min_real,
                "max_null_p90_fast": max_null,
                "strict_margin_fast": strict_margin,
                "strict_fast": strict_fast,
                "aux_score_fast": aux,
            }
        )
        if (i + 1) % 50 == 0:
            print(f"[QW-1985] fast search progress: {i + 1}/{total}", flush=True)

    fast_rows.sort(
        key=lambda x: (
            int(x["strict_fast"]),
            x["strict_margin_fast"],
            x["aux_score_fast"],
        ),
        reverse=True,
    )
    shortlist = fast_rows[:SHORTLIST_SIZE]
    print(f"[QW-1985] shortlist size: {len(shortlist)}", flush=True)

    final_rows = []
    for i, c in enumerate(shortlist):
        xi1 = float(c["xi1"])
        xi3 = float(c["xi3"])
        fold_full = []
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
                seed=155000 + i * 10 + f,
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
                "pass_count": pass_count,
                "min_real_full": min_real,
                "max_null_p90_full": max_null,
                "strict_margin_full": strict_margin,
                "aux_score_full": aux,
                "fold_results": fold_full,
            }
        )

    final_rows.sort(
        key=lambda x: (
            x["pass_count"],
            x["strict_margin_full"],
            x["aux_score_full"],
        ),
        reverse=True,
    )
    best = final_rows[0]

    verdict = (
        "STRICT_5OF5_PASS"
        if best["pass_count"] == 5
        else "STRICT_5OF5_NEAR_MISS"
        if best["pass_count"] == 4
        else "STRICT_5OF5_FAIL"
    )
    required = (
        "PROMOTE_TO_EXTERNAL_PREP_AND_LOCK_OPERATOR"
        if best["pass_count"] == 5
        else "ONE_FOLD_BORDERLINE_REQUIRES_OPERATOR_EXTENSION_OR_STRICTER_OBJECTIVE"
        if best["pass_count"] == 4
        else "EXTEND_OPERATOR_CLASS_AND_REPEAT"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1984.name, IN_QW1969.name],
        "search_config": {
            "grid_points": len(deltas) * len(deltas),
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
        "center": {"xi1": xi1_c, "xi3": xi3_c},
        "best": best,
        "top12": final_rows,
        "verdict": verdict,
        "required_next_step": required,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1985: STRICT 5/5 LOCAL PUSH",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Best Candidate",
        f"- xi1/xi3: {best['xi1']:.6f}/{best['xi3']:.6f}",
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

    print(f"[QW-1985] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1985] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1985] verdict={verdict} pass_count={best['pass_count']}/5")


if __name__ == "__main__":
    main()

