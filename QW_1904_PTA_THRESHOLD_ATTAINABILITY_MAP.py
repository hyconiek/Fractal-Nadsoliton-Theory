#!/usr/bin/env python3
"""
QW-1904: PTA threshold attainability map under feature-signal injection.

Purpose:
- Quantify how strong feature-linked signal must be to pass locked PTA V2 thresholds.
- Does NOT modify protocol; only stress-tests attainability on internal proxy data.
"""

from __future__ import annotations

import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1904_pta_threshold_attainability_map.json"
OUT_MD = ROOT / "RAPORT_QW1904_PTA_THRESHOLD_ATTAINABILITY_MAP.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    pta_path = ROOT / "external_confirmatory_v2" / "confirmatory_dataset_internal_proxy_wide" / "pta_v2_pairs.csv"
    if not pta_path.exists():
        raise RuntimeError(f"Missing PTA proxy dataset: {pta_path}")

    df0 = pd.read_csv(pta_path)

    d1850 = json.loads((ROOT / "report_qw1850_pta_v2_prereg_protocol.json").read_text(encoding="utf-8"))
    thr = d1850["protocol"]["pta_v2_protocol"]["thresholds"]

    mod1853 = load_module(ROOT / "QW_1853_JOINT_EXTERNAL_CONFIRMATORY_V2.py", "qw1853_mod_1904")

    # Feature score used only for stress scaling.
    z1 = (df0["f_autoc1"].to_numpy(dtype=float) - df0["f_autoc1"].mean()) / (df0["f_autoc1"].std() + 1e-12)
    z2 = (df0["f_switch"].to_numpy(dtype=float) - df0["f_switch"].mean()) / (df0["f_switch"].std() + 1e-12)
    z3 = (df0["f_std"].to_numpy(dtype=float) - df0["f_std"].mean()) / (df0["f_std"].std() + 1e-12)
    feature_score = 0.60 * z1 - 0.35 * z2 + 0.25 * z3

    # injected_hxy = clip(hxy + alpha * 0.05 * feature_score)
    # alpha=1 means 0.05-sigma feature contribution scale.
    alpha_grid = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0]

    rows = []
    first_pass_alpha = None

    for alpha in alpha_grid:
        df = df0.copy()
        df["hxy"] = np.clip(df["hxy"].to_numpy(dtype=float) + float(alpha) * 0.05 * feature_score, 0.0, 1.0)

        res = mod1853.eval_pta_v2(df, thresholds=thr)
        summ = res["summary"]
        flags = res["pass_flags"]
        all_pass = bool(all(bool(v) for v in flags.values()))

        row = {
            "alpha": float(alpha),
            "summary": {
                "mean_pair_mean_gain": float(summ["mean_pair_mean_gain"]),
                "bootstrap_lower95_mean_pair_mean_gain": float(summ["bootstrap_lower95_mean_pair_mean_gain"]),
                "prob_pair_mean_gain_positive": float(summ["prob_pair_mean_gain_positive"]),
                "one_sided_lower95_prob_pair_mean_gain_positive": float(summ["one_sided_lower95_prob_pair_mean_gain_positive"]),
            },
            "pass_flags": {k: bool(v) for k, v in flags.items()},
            "all_pass": all_pass,
        }
        rows.append(row)

        if all_pass and first_pass_alpha is None:
            first_pass_alpha = float(alpha)

    baseline = rows[0]

    if first_pass_alpha is not None and first_pass_alpha <= 2.0:
        verdict = "PTA_THRESHOLDS_ATTAINABLE_WITH_MODERATE_FEATURE_SIGNAL"
    elif first_pass_alpha is not None:
        verdict = "PTA_THRESHOLDS_ATTAINABLE_ONLY_WITH_STRONG_FEATURE_SIGNAL"
    else:
        verdict = "PTA_THRESHOLDS_NOT_ATTAINED_IN_TESTED_RANGE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_dataset": str(pta_path.relative_to(ROOT)),
        "thresholds": thr,
        "alpha_grid": alpha_grid,
        "results": rows,
        "baseline_alpha0": baseline,
        "first_pass_alpha": first_pass_alpha,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1904: PTA THRESHOLD ATTAINABILITY MAP",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- first_pass_alpha: {first_pass_alpha}",
        "",
        "## Baseline (alpha=0)",
        f"- mean_pair_mean_gain: {baseline['summary']['mean_pair_mean_gain']:.6f}",
        f"- bootstrap_lower95_mean_pair_mean_gain: {baseline['summary']['bootstrap_lower95_mean_pair_mean_gain']:.6f}",
        f"- prob_pair_mean_gain_positive: {baseline['summary']['prob_pair_mean_gain_positive']:.3f}",
        f"- one_sided_lower95_prob_pair_mean_gain_positive: {baseline['summary']['one_sided_lower95_prob_pair_mean_gain_positive']:.3f}",
        "",
        "## Thresholds",
        f"- mean_pair_mean_gain_min: {thr['mean_pair_mean_gain_min']}",
        f"- bootstrap_lower95_mean_pair_mean_gain_min: {thr['bootstrap_lower95_mean_pair_mean_gain_min']}",
        f"- prob_pair_mean_gain_positive_min: {thr['prob_pair_mean_gain_positive_min']}",
        f"- one_sided_lower95_prob_pair_mean_gain_positive_min: {thr['one_sided_lower95_prob_pair_mean_gain_positive_min']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1904] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1904] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
