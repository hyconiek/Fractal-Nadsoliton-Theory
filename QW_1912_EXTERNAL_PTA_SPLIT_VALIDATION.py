#!/usr/bin/env python3
"""
QW-1912: Discovery/holdout split validation for external PTA alpha family.

Purpose:
- reduce single-dataset selection bias by selecting alpha on discovery split only,
- verify locked PTA thresholds on disjoint holdout split.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1912_external_pta_split_validation.json"
OUT_MD = ROOT / "RAPORT_QW1912_EXTERNAL_PTA_SPLIT_VALIDATION.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_alpha_grid(text: str) -> List[float]:
    out = []
    for tok in text.split(","):
        s = tok.strip()
        if s:
            out.append(float(s))
    if not out:
        raise ValueError("alpha grid must not be empty")
    return out


def split_mask(pair_ids: pd.Series) -> np.ndarray:
    bits = []
    for pid in pair_ids.astype(str).tolist():
        h = hashlib.sha256(pid.encode("utf-8")).hexdigest()
        bits.append(int(h[-1], 16) % 2)
    return np.asarray(bits, dtype=int)


def feature_score(df: pd.DataFrame) -> np.ndarray:
    z1 = (df["f_autoc1"].to_numpy(dtype=float) - float(df["f_autoc1"].mean())) / (float(df["f_autoc1"].std()) + 1e-12)
    z2 = (df["f_switch"].to_numpy(dtype=float) - float(df["f_switch"].mean())) / (float(df["f_switch"].std()) + 1e-12)
    z3 = (df["f_std"].to_numpy(dtype=float) - float(df["f_std"].mean())) / (float(df["f_std"].std()) + 1e-12)
    return 0.60 * z1 - 0.35 * z2 + 0.25 * z3


def apply_alpha(df: pd.DataFrame, alpha: float, scale: float) -> pd.DataFrame:
    out = df.copy()
    score = feature_score(out)
    out["hxy"] = np.clip(out["hxy"].to_numpy(dtype=float) + float(alpha) * float(scale) * score, 0.0, 1.0)
    return out


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
        raise RuntimeError(f"Missing PTA source: {pta_path}")

    d1850 = json.loads((ROOT / "report_qw1850_pta_v2_prereg_protocol.json").read_text(encoding="utf-8"))
    thresholds = d1850["protocol"]["pta_v2_protocol"]["thresholds"]
    mod1853 = load_module(ROOT / "QW_1853_JOINT_EXTERNAL_CONFIRMATORY_V2.py", "qw1853_mod_1912")

    df = pd.read_csv(pta_path)
    if "pair_id" not in df.columns:
        raise RuntimeError("PTA file missing pair_id column required for split.")

    mask = split_mask(df["pair_id"])
    discovery = df[mask == 0].reset_index(drop=True)
    holdout = df[mask == 1].reset_index(drop=True)

    if len(discovery) < 200 or len(holdout) < 200:
        raise RuntimeError(f"Split too small: discovery={len(discovery)}, holdout={len(holdout)}")

    alphas = parse_alpha_grid(args.alpha_grid)

    discovery_rows = []
    selected_alpha = None

    for alpha in alphas:
        dfa = apply_alpha(discovery, alpha=alpha, scale=float(args.scale))
        res = mod1853.eval_pta_v2(dfa, thresholds=thresholds)
        flags = {k: bool(v) for k, v in res["pass_flags"].items()}
        all_pass = bool(all(flags.values()))
        discovery_rows.append(
            {
                "alpha": float(alpha),
                "summary": {k: float(v) for k, v in res["summary"].items()},
                "pass_flags": flags,
                "all_pass": all_pass,
            }
        )
        if all_pass and selected_alpha is None:
            selected_alpha = float(alpha)

    holdout_result = None
    if selected_alpha is not None:
        hfa = apply_alpha(holdout, alpha=selected_alpha, scale=float(args.scale))
        hres = mod1853.eval_pta_v2(hfa, thresholds=thresholds)
        holdout_result = {
            "alpha": selected_alpha,
            "summary": {k: float(v) for k, v in hres["summary"].items()},
            "pass_flags": {k: bool(v) for k, v in hres["pass_flags"].items()},
            "all_pass": bool(all(bool(v) for v in hres["pass_flags"].values())),
        }

    if selected_alpha is None:
        verdict = "SPLIT_VALIDATION_NO_ALPHA_PASS_ON_DISCOVERY"
    elif holdout_result is not None and holdout_result["all_pass"]:
        verdict = "SPLIT_VALIDATION_PASS_DISCOVERY_AND_HOLDOUT"
    else:
        verdict = "SPLIT_VALIDATION_FAIL_ON_HOLDOUT"

    out: Dict = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "pta_source": str(pta_path),
        "thresholds_locked": thresholds,
        "alpha_grid": [float(a) for a in alphas],
        "scale": float(args.scale),
        "split": {
            "method": "sha256(pair_id)_last_hex_parity",
            "n_discovery": int(len(discovery)),
            "n_holdout": int(len(holdout)),
        },
        "discovery_scan": discovery_rows,
        "selected_alpha_from_discovery": selected_alpha,
        "holdout_validation": holdout_result,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1912: EXTERNAL PTA SPLIT VALIDATION",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- PTA source: `{pta_path}`",
        f"- n_discovery: {out['split']['n_discovery']}",
        f"- n_holdout: {out['split']['n_holdout']}",
        f"- selected_alpha_from_discovery: {selected_alpha}",
        "",
        "## Holdout Result",
    ]

    if holdout_result is None:
        lines.append("- holdout not executed (no alpha passed discovery).")
    else:
        lines.extend(
            [
                f"- alpha: {holdout_result['alpha']}",
                f"- all_pass: {holdout_result['all_pass']}",
                f"- mean_pair_mean_gain: {holdout_result['summary']['mean_pair_mean_gain']:.6f}",
                f"- prob_pair_mean_gain_positive: {holdout_result['summary']['prob_pair_mean_gain_positive']:.3f}",
                f"- one_sided_lower95_prob_pair_mean_gain_positive: {holdout_result['summary']['one_sided_lower95_prob_pair_mean_gain_positive']:.3f}",
            ]
        )

    lines.extend(
        [
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1912] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1912] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1912] verdict={verdict}")


if __name__ == "__main__":
    main()

