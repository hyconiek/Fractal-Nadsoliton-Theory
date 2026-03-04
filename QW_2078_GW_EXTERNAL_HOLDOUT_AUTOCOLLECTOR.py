#!/usr/bin/env python3
"""
QW-2078: GW external holdout autocollector for QW-2077 validator.

Purpose:
- compute preregistered GW metrics from an external holdout feature table,
- fill observation JSON for QW-2077 without manual metric transcription,
- preserve source integrity metadata (sha256, row counts, pair counts).

This script does not assert "independent" by itself. Independence is a protocol property
of the data source and execution environment.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
PREREG_JSON = ROOT / "report_qw2076_empirical_prediction_preregistration.json"
TEMPLATE_JSON = ROOT / "empirical_observations_input_qw2077.template.json"
DEFAULT_INPUT = ROOT / "gw1831_window_features.csv"
DEFAULT_OBS_OUT = ROOT / "empirical_observations_input_qw2077.gw_autocollected.json"
OUT_JSON = ROOT / "report_qw2078_gw_external_holdout_autocollector.json"
OUT_MD = ROOT / "RAPORT_QW2078_GW_EXTERNAL_HOLDOUT_AUTOCOLLECTOR.md"

REQ_COLS = [
    "pair",
    "max_abs_corr",
    "mean_abs_corr",
    "corr_at_0ms",
    "corr_at_10ms",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rank_auc_pos_gt_neg(pos: np.ndarray, neg: np.ndarray) -> float:
    y = np.concatenate([np.ones(len(pos), dtype=int), np.zeros(len(neg), dtype=int)])
    s = np.concatenate([pos, neg])
    n1 = len(pos)
    n0 = len(neg)
    order = np.argsort(s)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    rs = float(np.sum(ranks[y == 1]))
    return float((rs - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def load_locked_gw_prediction() -> Dict:
    prereg = json.loads(PREREG_JSON.read_text(encoding="utf-8"))
    for pred in prereg["predictions"]:
        if pred.get("id") == "PRED-2076-GW-HOLDOUT":
            return pred
    raise KeyError("Missing PRED-2076-GW-HOLDOUT in preregistration JSON.")


def load_template() -> Dict:
    if TEMPLATE_JSON.exists():
        return json.loads(TEMPLATE_JSON.read_text(encoding="utf-8"))
    # Fallback minimal structure if template is absent.
    return {
        "pmns_cp": {
            "sin_delta_central": None,
            "sin_delta_ci95_low": None,
            "sin_delta_ci95_high": None,
            "source": "",
        },
        "cosmology_weff_nodes": [],
        "gw_future_holdout": {
            "auc_h1l1_vs_ctrl": None,
            "adv_shared_minus_ctrl_q90": None,
            "sep_median_h1l1_minus_ctrl": None,
            "control_median_gap": None,
            "n_windows": None,
            "source": "",
        },
        "notes": "",
    }


def main() -> None:
    data_path = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_INPUT.resolve()
    obs_out = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else DEFAULT_OBS_OUT.resolve()

    if not data_path.exists():
        raise FileNotFoundError(f"Input GW file not found: {data_path}")

    pred = load_locked_gw_prediction()
    locked_weights = pred["locked_weights"]
    thresholds = pred["hard_thresholds"]

    df = pd.read_csv(data_path)
    missing_cols = [c for c in REQ_COLS if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Input GW file missing required columns: {missing_cols}")

    pair = df["pair"].astype(str).to_numpy()
    score = (
        float(locked_weights["w_max_abs_corr"]) * df["max_abs_corr"].to_numpy(dtype=float)
        + float(locked_weights["w_mean_abs_corr"]) * df["mean_abs_corr"].to_numpy(dtype=float)
        + float(locked_weights["w_corr_at_0ms"]) * df["corr_at_0ms"].to_numpy(dtype=float)
        + float(locked_weights["w_corr_at_10ms"]) * df["corr_at_10ms"].to_numpy(dtype=float)
    )

    s_hl = score[pair == "H1-L1"]
    s_hv = score[pair == "H1-V1"]
    s_lv = score[pair == "L1-V1"]
    s_ctrl = np.concatenate([s_hv, s_lv])
    if len(s_hl) == 0 or len(s_hv) == 0 or len(s_lv) == 0:
        raise ValueError("Input GW file must contain pair labels: H1-L1, H1-V1, L1-V1.")

    q90 = float(np.quantile(s_ctrl, 0.90))
    metrics = {
        "auc_h1l1_vs_ctrl": float(rank_auc_pos_gt_neg(s_hl, s_ctrl)),
        "adv_shared_minus_ctrl_q90": float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90)),
        "sep_median_h1l1_minus_ctrl": float(np.median(s_hl) - np.median(s_ctrl)),
        "control_median_gap": float(abs(np.median(s_hv) - np.median(s_lv))),
    }
    checks = {
        "auc_ge_min": bool(metrics["auc_h1l1_vs_ctrl"] >= float(thresholds["auc_h1l1_vs_ctrl_min"])),
        "adv_ge_min": bool(metrics["adv_shared_minus_ctrl_q90"] >= float(thresholds["adv_shared_minus_ctrl_q90_min"])),
        "sep_ge_min": bool(metrics["sep_median_h1l1_minus_ctrl"] >= float(thresholds["sep_median_h1l1_minus_ctrl_min"])),
        "gap_le_max": bool(metrics["control_median_gap"] <= float(thresholds["control_median_gap_max"])),
    }

    obs = load_template()
    obs["gw_future_holdout"] = {
        **metrics,
        "n_windows": int(len(df)),
        "source": f"autocollected_from:{data_path.name}",
    }
    note = (
        f"QW-2078 autocollector used {data_path.name}; "
        f"sha256={sha256_file(data_path)}; rows={len(df)}."
    )
    if "gw1831_window_features.csv" in data_path.name:
        note += " WARNING: this is an internal legacy file, not a fresh independent holdout."
    obs["notes"] = ((obs.get("notes") or "").strip() + " " + note).strip()
    obs_out.write_text(json.dumps(obs, ensure_ascii=False, indent=2), encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": "GW_EXTERNAL_HOLDOUT_AUTOCOLLECTED",
        "input_file": str(data_path),
        "input_sha256": sha256_file(data_path),
        "input_rows": int(len(df)),
        "pair_counts": {k: int(v) for k, v in df["pair"].astype(str).value_counts().items()},
        "locked_weights": {k: float(v) for k, v in locked_weights.items()},
        "thresholds": {k: float(v) for k, v in thresholds.items()},
        "metrics": metrics,
        "checks": checks,
        "all_thresholds_pass": bool(all(checks.values())),
        "observation_output": str(obs_out),
        "required_next_step": (
            f"RUN_QW2077_WITH:{obs_out.name}"
        ),
        "independence_note": (
            "Independence must be ensured by external protocol (fresh holdout source + separate environment)."
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2078: GW EXTERNAL HOLDOUT AUTOCOLLECTOR",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Input: `{data_path}`",
        f"- SHA256: `{out['input_sha256']}`",
        f"- Rows: `{out['input_rows']}`",
        "",
        "## Metrics",
        f"- auc_h1l1_vs_ctrl: `{metrics['auc_h1l1_vs_ctrl']:.9f}`",
        f"- adv_shared_minus_ctrl_q90: `{metrics['adv_shared_minus_ctrl_q90']:.9f}`",
        f"- sep_median_h1l1_minus_ctrl: `{metrics['sep_median_h1l1_minus_ctrl']:.9f}`",
        f"- control_median_gap: `{metrics['control_median_gap']:.9f}`",
        "",
        "## Threshold Checks",
        f"- auc_ge_min: `{checks['auc_ge_min']}`",
        f"- adv_ge_min: `{checks['adv_ge_min']}`",
        f"- sep_ge_min: `{checks['sep_ge_min']}`",
        f"- gap_le_max: `{checks['gap_le_max']}`",
        f"- all_thresholds_pass: `{out['all_thresholds_pass']}`",
        "",
        "## Output",
        f"- observation JSON for QW-2077: `{obs_out}`",
        f"- report JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2078] Saved observation JSON: {obs_out}")
    print(f"[QW-2078] Saved report JSON: {OUT_JSON.name}")
    print(f"[QW-2078] Saved report MD:   {OUT_MD.name}")
    print(
        "[QW-2078] metrics="
        f"{metrics['auc_h1l1_vs_ctrl']:.4f}/"
        f"{metrics['adv_shared_minus_ctrl_q90']:.4f}/"
        f"{metrics['sep_median_h1l1_minus_ctrl']:.6f}/"
        f"{metrics['control_median_gap']:.6f}"
    )
    print(f"[QW-2078] all_thresholds_pass={out['all_thresholds_pass']}")


if __name__ == "__main__":
    main()

