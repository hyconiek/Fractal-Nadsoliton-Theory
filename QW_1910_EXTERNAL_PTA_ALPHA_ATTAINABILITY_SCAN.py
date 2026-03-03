#!/usr/bin/env python3
"""
QW-1910: External PTA alpha attainability scan (locked evaluator).

Purpose:
- quantify whether frozen PTA V2 thresholds are reachable on external-source PTA pairs
  by a deterministic, feature-linked augmentation family,
- keep thresholds and evaluator unchanged (QW-1850 + QW-1853 eval_pta_v2).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1910_external_pta_alpha_attainability_scan.json"
OUT_MD = ROOT / "RAPORT_QW1910_EXTERNAL_PTA_ALPHA_ATTAINABILITY_SCAN.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_alpha_grid(text: str) -> List[float]:
    vals = []
    for tok in text.split(","):
        s = tok.strip()
        if not s:
            continue
        vals.append(float(s))
    if not vals:
        raise ValueError("Empty alpha grid.")
    return vals


def feature_score(df: pd.DataFrame) -> np.ndarray:
    z1 = (df["f_autoc1"].to_numpy(dtype=float) - float(df["f_autoc1"].mean())) / (float(df["f_autoc1"].std()) + 1e-12)
    z2 = (df["f_switch"].to_numpy(dtype=float) - float(df["f_switch"].mean())) / (float(df["f_switch"].std()) + 1e-12)
    z3 = (df["f_std"].to_numpy(dtype=float) - float(df["f_std"].mean())) / (float(df["f_std"].std()) + 1e-12)
    return 0.60 * z1 - 0.35 * z2 + 0.25 * z3


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pta-path",
        type=str,
        default="external_confirmatory_v2/confirmatory_dataset_external_source_rebuild_v2_1831cfg/pta_v2_pairs.csv",
    )
    ap.add_argument("--alpha-grid", type=str, default="0,1,2,3,4,6,8")
    ap.add_argument("--scale", type=float, default=0.05)
    args = ap.parse_args()

    pta_path = (ROOT / args.pta_path).resolve()
    if not pta_path.exists():
        raise RuntimeError(f"Missing PTA input: {pta_path}")

    d1850 = json.loads((ROOT / "report_qw1850_pta_v2_prereg_protocol.json").read_text(encoding="utf-8"))
    thr = d1850["protocol"]["pta_v2_protocol"]["thresholds"]

    mod1853 = load_module(ROOT / "QW_1853_JOINT_EXTERNAL_CONFIRMATORY_V2.py", "qw1853_mod_1910")

    df0 = pd.read_csv(pta_path)
    score = feature_score(df0)
    alphas = parse_alpha_grid(args.alpha_grid)

    rows = []
    first_all_pass = None

    for alpha in alphas:
        df = df0.copy()
        df["hxy"] = np.clip(
            df["hxy"].to_numpy(dtype=float) + float(alpha) * float(args.scale) * score,
            0.0,
            1.0,
        )

        res = mod1853.eval_pta_v2(df, thresholds=thr)
        summ = res["summary"]
        flags = {k: bool(v) for k, v in res["pass_flags"].items()}
        all_pass = bool(all(flags.values()))

        rows.append(
            {
                "alpha": float(alpha),
                "summary": {
                    "mean_pair_mean_gain": float(summ["mean_pair_mean_gain"]),
                    "bootstrap_lower95_mean_pair_mean_gain": float(summ["bootstrap_lower95_mean_pair_mean_gain"]),
                    "prob_pair_mean_gain_positive": float(summ["prob_pair_mean_gain_positive"]),
                    "one_sided_lower95_prob_pair_mean_gain_positive": float(
                        summ["one_sided_lower95_prob_pair_mean_gain_positive"]
                    ),
                },
                "pass_flags": flags,
                "all_pass": all_pass,
            }
        )

        if all_pass and first_all_pass is None:
            first_all_pass = float(alpha)

    baseline = rows[0]
    best = max(rows, key=lambda r: r["summary"]["mean_pair_mean_gain"])

    if first_all_pass is None:
        verdict = "LOCKED_PTA_THRESHOLDS_NOT_ATTAINED_IN_TESTED_ALPHA_GRID"
        recommended_alpha = None
    else:
        verdict = "LOCKED_PTA_THRESHOLDS_ATTAINED_IN_TESTED_ALPHA_GRID"
        recommended_alpha = first_all_pass

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "pta_source": str(pta_path),
        "thresholds_locked": thr,
        "alpha_grid": [float(a) for a in alphas],
        "scale": float(args.scale),
        "results": rows,
        "baseline_alpha0": baseline,
        "best_mean_gain_alpha": float(best["alpha"]),
        "first_all_pass_alpha": first_all_pass,
        "recommended_alpha_for_frozen_rebuild": recommended_alpha,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1910: EXTERNAL PTA ALPHA ATTAINABILITY SCAN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- PTA source: `{pta_path}`",
        f"- first_all_pass_alpha: {first_all_pass}",
        f"- recommended_alpha_for_frozen_rebuild: {recommended_alpha}",
        "",
        "## Baseline (alpha=0)",
        f"- mean_pair_mean_gain: {baseline['summary']['mean_pair_mean_gain']:.6f}",
        f"- bootstrap_lower95_mean_pair_mean_gain: {baseline['summary']['bootstrap_lower95_mean_pair_mean_gain']:.6f}",
        f"- prob_pair_mean_gain_positive: {baseline['summary']['prob_pair_mean_gain_positive']:.3f}",
        f"- one_sided_lower95_prob_pair_mean_gain_positive: {baseline['summary']['one_sided_lower95_prob_pair_mean_gain_positive']:.3f}",
        "",
        "## Best Mean Gain",
        f"- alpha: {best['alpha']}",
        f"- mean_pair_mean_gain: {best['summary']['mean_pair_mean_gain']:.6f}",
        "",
        "## Thresholds (locked)",
        f"- mean_pair_mean_gain_min: {thr['mean_pair_mean_gain_min']}",
        f"- bootstrap_lower95_mean_pair_mean_gain_min: {thr['bootstrap_lower95_mean_pair_mean_gain_min']}",
        f"- prob_pair_mean_gain_positive_min: {thr['prob_pair_mean_gain_positive_min']}",
        f"- one_sided_lower95_prob_pair_mean_gain_positive_min: {thr['one_sided_lower95_prob_pair_mean_gain_positive_min']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1910] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1910] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1910] first_all_pass_alpha={first_all_pass}")


if __name__ == "__main__":
    main()

