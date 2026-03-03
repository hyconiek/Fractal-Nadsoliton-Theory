#!/usr/bin/env python3
"""
QW-1831: Event-windowed GW coherent feature extractor.

Builds detector-pair features from local unfiltered strain files (H1/L1/V1)
for GPS 1266965117.

Outputs:
- per-window feature table (CSV),
- aggregated per-pair diagnostics (JSON/MD).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, detrend


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "raw_strain_unfiltered"
OUT_CSV = ROOT / "gw1831_window_features.csv"
OUT_JSON = ROOT / "report_qw1831_event_windowed_gw_features.json"
OUT_MD = ROOT / "RAPORT_QW1831_EVENT_WINDOWED_GW_FEATURES.md"

FS = 4096
GPS = 1266965117
WINDOW_SEC = 8.0
STEP_SEC = 2.0
LAG_MAX_MS = 15.0
BAND = (20.0, 400.0)


def load_strain(det: str) -> np.ndarray:
    fp = DATA_DIR / f"{det}_unfiltered_{GPS}.h5"
    if not fp.exists():
        raise FileNotFoundError(f"Missing file: {fp}")
    with h5py.File(fp, "r") as f:
        x = np.asarray(f["strain"][:], dtype=float)
    return x


def bandpass(x: np.ndarray, fs: int, f_lo: float, f_hi: float, order: int = 4) -> np.ndarray:
    nyq = 0.5 * fs
    b, a = butter(order, [f_lo / nyq, f_hi / nyq], btype="band")
    return filtfilt(b, a, x)


def zscore(x: np.ndarray) -> np.ndarray:
    s = float(np.std(x))
    # GW strain amplitude can be ~1e-21, so use near-zero test instead of fixed macroscopic threshold.
    if (not np.isfinite(s)) or s <= 1e-30:
        return np.zeros_like(x)
    return (x - float(np.mean(x))) / s


def pair_window_features(x: np.ndarray, y: np.ndarray, fs: int, lag_max_ms: float) -> Dict[str, float]:
    xw = zscore(x)
    yw = zscore(y)

    max_lag = int(round((lag_max_ms / 1000.0) * fs))
    lags = np.arange(-max_lag, max_lag + 1)

    corrs = []
    for lag in lags:
        if lag < 0:
            a = xw[-lag:]
            b = yw[: len(yw) + lag]
        elif lag > 0:
            a = xw[: len(xw) - lag]
            b = yw[lag:]
        else:
            a = xw
            b = yw
        if len(a) < 10:
            corrs.append(np.nan)
            continue
        corrs.append(float(np.corrcoef(a, b)[0, 1]))

    c = np.array(corrs, dtype=float)
    if np.all(~np.isfinite(c)):
        return {
            "max_abs_corr": np.nan,
            "best_lag_ms": np.nan,
            "corr_at_0ms": np.nan,
            "corr_at_10ms": np.nan,
            "mean_abs_corr": np.nan,
        }

    idx = int(np.nanargmax(np.abs(c)))
    best_lag_ms = (lags[idx] / fs) * 1000.0

    lag10 = int(round(0.010 * fs))
    idx10 = np.where(lags == lag10)[0]
    c10 = float(c[idx10[0]]) if len(idx10) else np.nan

    idx0 = np.where(lags == 0)[0]
    c0 = float(c[idx0[0]]) if len(idx0) else np.nan

    return {
        "max_abs_corr": float(np.nanmax(np.abs(c))),
        "best_lag_ms": float(best_lag_ms),
        "corr_at_0ms": c0,
        "corr_at_10ms": c10,
        "mean_abs_corr": float(np.nanmean(np.abs(c))),
    }


def summarize_pair(df: pd.DataFrame) -> Dict[str, float]:
    out = {
        "n_windows": int(len(df)),
    }
    for col in ["max_abs_corr", "best_lag_ms", "corr_at_0ms", "corr_at_10ms", "mean_abs_corr"]:
        v = df[col].to_numpy(dtype=float)
        out[f"{col}_mean"] = float(np.nanmean(v))
        out[f"{col}_median"] = float(np.nanmedian(v))
        out[f"{col}_std"] = float(np.nanstd(v))
        out[f"{col}_p10"] = float(np.nanquantile(v, 0.10))
        out[f"{col}_p90"] = float(np.nanquantile(v, 0.90))

    out["prob_corr10_gt_0p02"] = float(np.mean(df["corr_at_10ms"].to_numpy(dtype=float) > 0.02))
    out["prob_maxabs_gt_0p02"] = float(np.mean(df["max_abs_corr"].to_numpy(dtype=float) > 0.02))
    out["prob_bestlag_abs_lt_3ms"] = float(np.mean(np.abs(df["best_lag_ms"].to_numpy(dtype=float)) < 3.0))
    return out


def main() -> None:
    x_h1 = load_strain("H1")
    x_l1 = load_strain("L1")
    x_v1 = load_strain("V1")

    # harmonize length
    n = min(len(x_h1), len(x_l1), len(x_v1))
    x_h1 = x_h1[:n]
    x_l1 = x_l1[:n]
    x_v1 = x_v1[:n]

    # preprocessing
    def prep(x: np.ndarray) -> np.ndarray:
        y = detrend(x)
        y = bandpass(y, fs=FS, f_lo=BAND[0], f_hi=BAND[1], order=4)
        return zscore(y)

    p_h1 = prep(x_h1)
    p_l1 = prep(x_l1)
    p_v1 = prep(x_v1)

    pairs: List[Tuple[str, np.ndarray, np.ndarray]] = [
        ("H1-L1", p_h1, p_l1),
        ("H1-V1", p_h1, p_v1),
        ("L1-V1", p_l1, p_v1),
    ]

    w = int(round(WINDOW_SEC * FS))
    step = int(round(STEP_SEC * FS))
    rows: List[Dict[str, float]] = []

    for pair_name, a, b in pairs:
        k = 0
        for start in range(0, n - w + 1, step):
            aa = a[start : start + w]
            bb = b[start : start + w]
            feats = pair_window_features(aa, bb, fs=FS, lag_max_ms=LAG_MAX_MS)
            row = {
                "pair": pair_name,
                "window_idx": k,
                "start": int(start),
                "stop": int(start + w),
            }
            row.update(feats)
            rows.append(row)
            k += 1

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)

    pair_summary = {}
    for pair_name in sorted(df["pair"].unique()):
        pair_summary[pair_name] = summarize_pair(df[df["pair"] == pair_name].copy())

    # cross-pair consistency diagnostics
    med_corr10 = [pair_summary[p]["corr_at_10ms_median"] for p in pair_summary]
    med_maxabs = [pair_summary[p]["max_abs_corr_median"] for p in pair_summary]

    consistency = {
        "std_pair_median_corr10": float(np.std(med_corr10)),
        "std_pair_median_maxabs": float(np.std(med_maxabs)),
        "range_pair_median_corr10": float(np.max(med_corr10) - np.min(med_corr10)),
        "range_pair_median_maxabs": float(np.max(med_maxabs) - np.min(med_maxabs)),
    }

    # Prototype pass flags for GW redesign stage (not final physics claim gate)
    pass_windows = all(pair_summary[p]["n_windows"] >= 100 for p in pair_summary)
    pass_lag_signal = pair_summary["H1-L1"]["prob_corr10_gt_0p02"] >= 0.10
    pass_pair_consistency = consistency["std_pair_median_maxabs"] <= 0.01

    if pass_windows and pass_lag_signal and pass_pair_consistency:
        verdict = "EVENT_WINDOW_FEATURE_BASELINE_SUPPORTED"
    elif pass_windows:
        verdict = "EVENT_WINDOW_FEATURE_BASELINE_PARTIAL"
    else:
        verdict = "EVENT_WINDOW_FEATURE_BASELINE_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "config": {
            "gps": GPS,
            "fs": FS,
            "window_sec": WINDOW_SEC,
            "step_sec": STEP_SEC,
            "lag_max_ms": LAG_MAX_MS,
            "band_hz": [BAND[0], BAND[1]],
            "n_samples": int(n),
        },
        "pair_summary": pair_summary,
        "consistency": consistency,
        "pass_flags": {
            "enough_windows": bool(pass_windows),
            "lag_signal_presence_h1l1": bool(pass_lag_signal),
            "inter_pair_consistency": bool(pass_pair_consistency),
        },
        "verdict": verdict,
        "artifacts": {
            "csv": OUT_CSV.name,
        },
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1831: EVENT-WINDOWED GW FEATURES",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- GPS: {GPS}",
        f"- n samples: {n}",
        f"- Window/step [s]: {WINDOW_SEC}/{STEP_SEC}",
        f"- Band [Hz]: {BAND[0]}-{BAND[1]}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pair Summary (selected)",
    ]
    for p in ["H1-L1", "H1-V1", "L1-V1"]:
        s = pair_summary[p]
        lines.append(
            f"- {p}: n={s['n_windows']}, median corr@10ms={s['corr_at_10ms_median']:.6f}, "
            f"median max|corr|={s['max_abs_corr_median']:.6f}, P(corr10>0.02)={s['prob_corr10_gt_0p02']:.3f}"
        )

    lines.extend(
        [
            "",
            "## Consistency",
            f"- std pair median corr@10ms: {consistency['std_pair_median_corr10']:.6f}",
            f"- std pair median max|corr|: {consistency['std_pair_median_maxabs']:.6f}",
            "",
            "## Pass Flags",
            f"- enough_windows: {pass_windows}",
            f"- lag_signal_presence_h1l1: {pass_lag_signal}",
            f"- inter_pair_consistency: {pass_pair_consistency}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
            f"- CSV: `{OUT_CSV.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1831] Saved CSV:  {OUT_CSV.name}")
    print(f"[QW-1831] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1831] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
