#!/usr/bin/env python3
"""
QW-1991: Fold-2 random-null suppression search (frozen tri-basis).
Primary target: reduce fold-2 random-null p90 while preserving global hard-lock metrics.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import ROOT
from QW_1989_CONSTRAINED_ADVERSARIAL_AUDIT import constrained_adversarial_rate
from QW_1990_HARD_LOCK_RANDOM_NULL_SEARCH import (
    TARGET_NULL_RANDOM_P90_MAX,
    TARGET_REAL_BLOCK_MIN,
    TARGET_REAL_IID_MIN,
    eval_candidate,
    unique_rows,
)
from QW_1986_TRI_BASIS_STRICT_5OF5_ATTEMPT import build_fold_channels_tri
from QW_1988_TRI_BASIS_HARD_LOCK_STRESS import bootstrap_pass_rate_block


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1990 = ROOT / "report_qw1990_hard_lock_random_null_search.json"
OUT_JSON = ROOT / "report_qw1991_fold2_random_null_suppression_search.json"
OUT_MD = ROOT / "RAPORT_QW1991_FOLD2_RANDOM_NULL_SUPPRESSION_SEARCH.md"

FAST_REAL_IID_BOOT = 130
FAST_NULL_TRIALS = 6
FAST_NULL_BOOT = 30

FULL_REAL_IID_BOOT = 1400
FULL_REAL_BLOCK_BOOT = 900
FULL_REAL_BLOCK_LEN = 10
FULL_NULL_TRIALS = 18
FULL_NULL_BOOT = 80
FULL_ADV_FLIP_BUDGET = 4

SHORTLIST_SIZE = 10


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    r1990 = json.loads(IN_QW1990.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]

    center = {
        "xi1": float(r1990["best"]["xi1"]),
        "xi3": float(r1990["best"]["xi3"]),
        "xi4": float(r1990["best"]["xi4"]),
    }

    d13 = np.array([-0.00010, -0.00005, 0.0, 0.00005, 0.00010], dtype=float)
    d4 = np.array([-0.00016, -0.00008, 0.0, 0.00008, 0.00016], dtype=float)
    rows = []
    for a in d13:
        for b in d13:
            for c in d4:
                rows.append(
                    {
                        "xi1": float(np.clip(center["xi1"] + a, 0.0001, 0.0032)),
                        "xi3": float(np.clip(center["xi3"] + b, -0.0002, 0.0032)),
                        "xi4": float(np.clip(center["xi4"] + c, -0.0032, 0.0032)),
                    }
                )
    rng = np.random.default_rng(210000)
    for _ in range(120):
        rows.append(
            {
                "xi1": float(np.clip(rng.normal(center["xi1"], 0.00009), 0.0001, 0.0032)),
                "xi3": float(np.clip(rng.normal(center["xi3"], 0.00009), -0.0002, 0.0032)),
                "xi4": float(np.clip(rng.normal(center["xi4"], 0.00012), -0.0032, 0.0032)),
            }
        )
    candidates = unique_rows(rows)

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
        fr = []
        for f, dff in enumerate(fold_dfs):
            row = eval_candidate(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                real_iid_boot=FAST_REAL_IID_BOOT,
                null_trials=FAST_NULL_TRIALS,
                null_boot=FAST_NULL_BOOT,
                seed=211000 + i * 10 + f,
            )
            fr.append({"fold": f, **row})

        min_real = float(min(r["real_iid"] for r in fr))
        max_null = float(max(r["null_random_p90"] for r in fr))
        fold2_null = float([r["null_random_p90"] for r in fr if r["fold"] == 2][0])
        hard_margin = float(min(min_real - TARGET_REAL_IID_MIN, TARGET_NULL_RANDOM_P90_MAX - max_null))
        # Targeted objective: prioritize fold-2 suppression with global safety.
        targeted = float((TARGET_NULL_RANDOM_P90_MAX - fold2_null) + 0.70 * hard_margin)
        aux = float(np.mean([r["real_iid"] for r in fr]) - np.mean([r["null_random_p90"] for r in fr]))

        fast_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "xi4": xi4,
                "fold_fast": fr,
                "min_real_iid_fast": min_real,
                "max_null_random_p90_fast": max_null,
                "fold2_null_random_p90_fast": fold2_null,
                "hard_margin_fast": hard_margin,
                "targeted_score_fast": targeted,
                "aux_fast": aux,
            }
        )
        if (i + 1) % 40 == 0:
            print(f"[QW-1991] fast search progress: {i + 1}/{total}", flush=True)

    fast_rows.sort(
        key=lambda x: (
            x["targeted_score_fast"],
            x["hard_margin_fast"],
            x["aux_fast"],
        ),
        reverse=True,
    )
    shortlist = fast_rows[:SHORTLIST_SIZE]
    print(f"[QW-1991] shortlist size: {len(shortlist)}", flush=True)

    final_rows = []
    for i, c in enumerate(shortlist):
        xi1 = float(c["xi1"])
        xi3 = float(c["xi3"])
        xi4 = float(c["xi4"])

        fr = []
        for f, dff in enumerate(fold_dfs):
            base = eval_candidate(
                dff=dff,
                kernel=kernel,
                params=params,
                thr=thr,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                real_iid_boot=FULL_REAL_IID_BOOT,
                null_trials=FULL_NULL_TRIALS,
                null_boot=FULL_NULL_BOOT,
                seed=215000 + i * 10 + f,
            )

            s_hl, s_hv, s_lv, pairs, c1, c3, c4 = build_fold_channels_tri(dff, kernel, params, xi1, xi3, xi4)
            real_block = bootstrap_pass_rate_block(
                s_hl=s_hl,
                s_hv=s_hv,
                s_lv=s_lv,
                thr=thr,
                n_boot=FULL_REAL_BLOCK_BOOT,
                seed=216000 + i * 10 + f,
                block_len=FULL_REAL_BLOCK_LEN,
            )

            pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
            pairs_vec = dff["pair"].map(pair_map).to_numpy(dtype=int)
            s_full = np.zeros(len(pairs_vec), dtype=float)
            s_full[pairs_vec == 0] = s_hl
            s_full[pairs_vec == 1] = s_hv
            s_full[pairs_vec == 2] = s_lv
            s_base = s_full - xi1 * c1 - xi3 * c3 - xi4 * c4
            ctrl_idx = np.where(pairs_vec != 0)[0]
            control_order = np.argsort(dff.iloc[ctrl_idx]["window_idx"].to_numpy(dtype=int))
            adv_rate = constrained_adversarial_rate(
                s_base=s_base,
                pairs=pairs_vec,
                c1=c1,
                c3=c3,
                c4=c4,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                thr=thr,
                control_order=control_order,
                max_flips=FULL_ADV_FLIP_BUDGET,
                seed=217000 + i * 10 + f,
            )

            fr.append(
                {
                    "fold": f,
                    "real_iid_full": base["real_iid"],
                    "real_block_full": real_block,
                    "null_random_mean_full": base["null_random_mean"],
                    "null_random_p90_full": base["null_random_p90"],
                    "adv_constrained_full": adv_rate,
                }
            )

        min_real_iid = float(min(r["real_iid_full"] for r in fr))
        min_real_block = float(min(r["real_block_full"] for r in fr))
        max_null_p90 = float(max(r["null_random_p90_full"] for r in fr))
        fold2_null_p90 = float([r["null_random_p90_full"] for r in fr if r["fold"] == 2][0])
        max_adv = float(max(r["adv_constrained_full"] for r in fr))
        hard_margin = float(min(min_real_iid - TARGET_REAL_IID_MIN, TARGET_NULL_RANDOM_P90_MAX - max_null_p90))
        targeted = float((TARGET_NULL_RANDOM_P90_MAX - fold2_null_p90) + 0.70 * hard_margin)
        hard_pass = bool(
            min_real_iid >= TARGET_REAL_IID_MIN
            and min_real_block >= TARGET_REAL_BLOCK_MIN
            and max_null_p90 <= TARGET_NULL_RANDOM_P90_MAX
            and max_adv <= 0.45
        )
        aux = float(np.mean([r["real_iid_full"] for r in fr]) - np.mean([r["null_random_p90_full"] for r in fr]))

        final_rows.append(
            {
                "xi1": xi1,
                "xi3": xi3,
                "xi4": xi4,
                "hard_pass": hard_pass,
                "min_real_iid_full": min_real_iid,
                "min_real_block_full": min_real_block,
                "max_null_random_p90_full": max_null_p90,
                "fold2_null_random_p90_full": fold2_null_p90,
                "max_adv_constrained_full": max_adv,
                "hard_margin_full": hard_margin,
                "targeted_score_full": targeted,
                "aux_full": aux,
                "fold_results": fr,
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
    verdict = "FOLD2_SUPPRESSION_HARD_PASS" if best["hard_pass"] else "FOLD2_SUPPRESSION_HARD_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1990.name],
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
        "# RAPORT QW-1991: FOLD-2 RANDOM-NULL SUPPRESSION SEARCH",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Best Candidate",
        f"- xi1/xi3/xi4: {best['xi1']:.6f}/{best['xi3']:.6f}/{best['xi4']:.6f}",
        f"- hard_pass: {best['hard_pass']}",
        f"- min_real_iid_full: {100.0 * best['min_real_iid_full']:.2f}%",
        f"- min_real_block_full: {100.0 * best['min_real_block_full']:.2f}%",
        f"- max_null_random_p90_full: {100.0 * best['max_null_random_p90_full']:.2f}%",
        f"- fold2_null_random_p90_full: {100.0 * best['fold2_null_random_p90_full']:.2f}%",
        f"- max_adv_constrained_full: {100.0 * best['max_adv_constrained_full']:.2f}%",
        f"- hard_margin_full: {100.0 * best['hard_margin_full']:.2f} pp",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1991] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1991] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1991] verdict={verdict}")


if __name__ == "__main__":
    main()
