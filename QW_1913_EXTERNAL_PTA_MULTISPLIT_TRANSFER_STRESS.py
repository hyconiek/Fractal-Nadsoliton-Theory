#!/usr/bin/env python3
"""
QW-1913: External PTA multisplit transfer stress test.

Goal:
- test whether alpha selected on discovery partitions transfers to disjoint holdouts
  across multiple splits,
- keep evaluator and thresholds frozen (QW-1853 eval_pta_v2 + QW-1850 thresholds).
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
OUT_JSON = ROOT / "report_qw1913_external_pta_multisplit_transfer_stress.json"
OUT_MD = ROOT / "RAPORT_QW1913_EXTERNAL_PTA_MULTISPLIT_TRANSFER_STRESS.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_alpha_grid(text: str) -> List[float]:
    vals: List[float] = []
    for tok in text.split(","):
        s = tok.strip()
        if s:
            vals.append(float(s))
    if not vals:
        raise ValueError("Empty alpha grid.")
    return vals


def split_index(pair_id: str, k_folds: int) -> int:
    h = hashlib.sha256(pair_id.encode("utf-8")).hexdigest()
    return int(h[-8:], 16) % k_folds


def feature_score(df: pd.DataFrame) -> np.ndarray:
    z1 = (df["f_autoc1"].to_numpy(dtype=float) - float(df["f_autoc1"].mean())) / (float(df["f_autoc1"].std()) + 1e-12)
    z2 = (df["f_switch"].to_numpy(dtype=float) - float(df["f_switch"].mean())) / (float(df["f_switch"].std()) + 1e-12)
    z3 = (df["f_std"].to_numpy(dtype=float) - float(df["f_std"].mean())) / (float(df["f_std"].std()) + 1e-12)
    return 0.60 * z1 - 0.35 * z2 + 0.25 * z3


def apply_alpha(df: pd.DataFrame, alpha: float, scale: float) -> pd.DataFrame:
    out = df.copy()
    s = feature_score(out)
    out["hxy"] = np.clip(out["hxy"].to_numpy(dtype=float) + float(alpha) * float(scale) * s, 0.0, 1.0)
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
    ap.add_argument("--k-folds", type=int, default=3)
    args = ap.parse_args()

    if args.k_folds < 2:
        raise ValueError("k-folds must be >= 2.")

    pta_path = (ROOT / args.pta_path).resolve()
    if not pta_path.exists():
        raise RuntimeError(f"Missing PTA source: {pta_path}")

    d1850 = json.loads((ROOT / "report_qw1850_pta_v2_prereg_protocol.json").read_text(encoding="utf-8"))
    thresholds = d1850["protocol"]["pta_v2_protocol"]["thresholds"]
    mod1853 = load_module(ROOT / "QW_1853_JOINT_EXTERNAL_CONFIRMATORY_V2.py", "qw1853_mod_1913")

    df = pd.read_csv(pta_path)
    if "pair_id" not in df.columns:
        raise RuntimeError("PTA source requires pair_id column.")

    fold_ids = np.array([split_index(str(x), args.k_folds) for x in df["pair_id"]], dtype=int)
    alphas = parse_alpha_grid(args.alpha_grid)

    fold_rows: List[Dict] = []
    holdout_passes = 0
    folds_with_selected = 0

    for fold in range(args.k_folds):
        discovery = df[fold_ids != fold].reset_index(drop=True)
        holdout = df[fold_ids == fold].reset_index(drop=True)

        if len(discovery) < 400 or len(holdout) < 200:
            fold_rows.append(
                {
                    "fold": int(fold),
                    "n_discovery": int(len(discovery)),
                    "n_holdout": int(len(holdout)),
                    "status": "SKIPPED_SPLIT_TOO_SMALL",
                }
            )
            continue

        scan_rows = []
        selected_alpha = None

        for alpha in alphas:
            dfa = apply_alpha(discovery, alpha=alpha, scale=float(args.scale))
            res = mod1853.eval_pta_v2(dfa, thresholds=thresholds)
            flags = {k: bool(v) for k, v in res["pass_flags"].items()}
            all_pass = bool(all(flags.values()))
            scan_rows.append(
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
            folds_with_selected += 1
            hfa = apply_alpha(holdout, alpha=selected_alpha, scale=float(args.scale))
            hres = mod1853.eval_pta_v2(hfa, thresholds=thresholds)
            holdout_result = {
                "alpha": selected_alpha,
                "summary": {k: float(v) for k, v in hres["summary"].items()},
                "pass_flags": {k: bool(v) for k, v in hres["pass_flags"].items()},
                "all_pass": bool(all(bool(v) for v in hres["pass_flags"].values())),
            }
            if holdout_result["all_pass"]:
                holdout_passes += 1

        fold_rows.append(
            {
                "fold": int(fold),
                "n_discovery": int(len(discovery)),
                "n_holdout": int(len(holdout)),
                "scan": scan_rows,
                "selected_alpha": selected_alpha,
                "holdout": holdout_result,
                "status": "OK",
            }
        )

    selected_alphas = [float(r["selected_alpha"]) for r in fold_rows if r.get("selected_alpha") is not None]
    if folds_with_selected > 0:
        holdout_pass_rate = float(holdout_passes / folds_with_selected)
    else:
        holdout_pass_rate = 0.0

    if folds_with_selected == args.k_folds and holdout_pass_rate >= 1.0:
        verdict = "MULTISPLIT_TRANSFER_PASS_ALL_FOLDS"
    elif folds_with_selected >= max(1, args.k_folds - 1) and holdout_pass_rate >= 0.67:
        verdict = "MULTISPLIT_TRANSFER_PARTIAL"
    else:
        verdict = "MULTISPLIT_TRANSFER_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "pta_source": str(pta_path),
        "thresholds_locked": thresholds,
        "alpha_grid": [float(a) for a in alphas],
        "scale": float(args.scale),
        "k_folds": int(args.k_folds),
        "fold_results": fold_rows,
        "summary": {
            "folds_with_selected_alpha": int(folds_with_selected),
            "holdout_passes": int(holdout_passes),
            "holdout_pass_rate_given_selected": holdout_pass_rate,
            "selected_alphas": selected_alphas,
            "selected_alpha_median": float(np.median(selected_alphas)) if selected_alphas else None,
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1913: EXTERNAL PTA MULTISPLIT TRANSFER STRESS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- k_folds: {args.k_folds}",
        f"- folds_with_selected_alpha: {out['summary']['folds_with_selected_alpha']}",
        f"- holdout_passes: {out['summary']['holdout_passes']}",
        f"- holdout_pass_rate_given_selected: {out['summary']['holdout_pass_rate_given_selected']:.3f}",
        f"- selected_alphas: {out['summary']['selected_alphas']}",
        f"- selected_alpha_median: {out['summary']['selected_alpha_median']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1913] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1913] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1913] verdict={verdict}")


if __name__ == "__main__":
    main()

