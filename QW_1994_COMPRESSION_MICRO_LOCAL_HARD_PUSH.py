#!/usr/bin/env python3
"""
QW-1994: Micro-local hard push around QW-1993 compression candidate.
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
IN_QW1993 = ROOT / "report_qw1993_global_score_compression_hard_lock_search.json"
OUT_JSON = ROOT / "report_qw1994_compression_micro_local_hard_push.json"
OUT_MD = ROOT / "RAPORT_QW1994_COMPRESSION_MICRO_LOCAL_HARD_PUSH.md"

FAST_REAL_IID_BOOT = 180
FAST_NULL_TRIALS = 8
FAST_NULL_BOOT = 40

FULL_REAL_IID_BOOT = 1800
FULL_REAL_BLOCK_BOOT = 1100
FULL_REAL_BLOCK_LEN = 10
FULL_NULL_TRIALS = 24
FULL_NULL_BOOT = 100
FULL_ADV_FLIP_BUDGET = 4

SHORTLIST_SIZE = 12


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    r1993 = json.loads(IN_QW1993.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]
    center = r1993["best"]

    rows = []
    d13 = np.array([-0.00002, -0.00001, 0.0, 0.00001, 0.00002], dtype=float)
    d4 = np.array([-0.00003, 0.0, 0.00003], dtype=float)
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

    rng = np.random.default_rng(240000)
    for _ in range(140):
        rows.append(
            {
                "xi1": float(np.clip(rng.normal(float(center["xi1"]), 0.000015), 0.0001, 0.0032)),
                "xi3": float(np.clip(rng.normal(float(center["xi3"]), 0.000015), -0.0002, 0.0032)),
                "xi4": float(np.clip(rng.normal(float(center["xi4"]), 0.000025), -0.0032, 0.0032)),
                "gamma_c": float(np.clip(rng.normal(float(center["gamma_c"]), 0.010), 0.78, 1.02)),
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
                seed=241000 + i * 10 + f,
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
        if (i + 1) % 80 == 0:
            print(f"[QW-1994] fast search progress: {i + 1}/{total}", flush=True)

    fast_rows.sort(key=lambda x: (x["hard_margin_fast"], x["aux_fast"]), reverse=True)
    shortlist = fast_rows[:SHORTLIST_SIZE]
    print(f"[QW-1994] shortlist size: {len(shortlist)}", flush=True)

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
                seed=245000 + i * 10 + f,
            )
            # for block metric reuse compressed real channels through helper internals
            # we estimate block robustness by re-running eval at block-level using constructed scores indirectly:
            # evaluate with synthetic bootstrap on the same compressed score distribution
            # (consistent with QW-1993 methodology).
            from QW_1993_GLOBAL_SCORE_COMPRESSION_HARD_LOCK_SEARCH import build_components, compress_score

            base_score, c1, c3, c4, pairs = build_components(dff, kernel, params)
            s_raw = base_score + xi1 * c1 + xi3 * c3 + xi4 * c4
            s = compress_score(s_raw, gc)
            s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]
            real_block = bootstrap_pass_rate_block(
                s_hl=s_hl,
                s_hv=s_hv,
                s_lv=s_lv,
                thr=thr,
                n_boot=FULL_REAL_BLOCK_BOOT,
                seed=246000 + i * 10 + f,
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
                seed=247000 + i * 10 + f,
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
    verdict = "COMPRESSION_MICRO_HARD_PASS" if best["hard_pass"] else "COMPRESSION_MICRO_HARD_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1993.name],
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
        },
        "best": best,
        "top12": final_rows,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1994: COMPRESSION MICRO-LOCAL HARD PUSH",
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

    print(f"[QW-1994] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1994] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1994] verdict={verdict}")


if __name__ == "__main__":
    main()

