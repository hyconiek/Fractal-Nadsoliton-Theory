#!/usr/bin/env python3
"""
Phase 65 strict repair:
- canonical whitening + robust lag profile,
- permutation significance for max-abs lag correlation.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from qw1660_strict_common import exact_phase_only_whiten, read_strain, safe_pearson


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw_strain_unfiltered"
OUT = ROOT / "QW_1660_v65_MicroTimeShift_strict.json"


def lag_profile(x: np.ndarray, y: np.ndarray, fs: int, max_ms: float = 25.0) -> dict:
    max_lag = int(round(max_ms * fs / 1000.0))
    lags = np.arange(-max_lag, max_lag + 1)
    center_start = max_lag
    center_end = len(y) - max_lag
    y0 = y[center_start:center_end]
    corr = []
    for lag in lags:
        xs = x[center_start + lag : center_end + lag]
        corr.append(safe_pearson(xs, y0))
    corr = np.array(corr, dtype=float)
    idx_abs = int(np.argmax(np.abs(corr)))
    idx_pos = int(np.argmax(corr))
    return {
        "lags_samples": lags,
        "corr": corr,
        "best_abs_lag_samples": int(lags[idx_abs]),
        "best_abs_corr": float(corr[idx_abs]),
        "best_pos_lag_samples": int(lags[idx_pos]),
        "best_pos_corr": float(corr[idx_pos]),
    }


def main() -> None:
    gps = 1266965117
    fs = 4096
    n_samples = int(os.getenv("PHASE65_STRICT_N_SAMPLES", "524288"))
    n_perm = int(os.getenv("PHASE65_STRICT_N_PERM", "200"))
    seed = int(os.getenv("PHASE65_STRICT_SEED", "6501"))
    rng = np.random.default_rng(seed)

    x_h1 = read_strain(RAW, "H1", gps, n_samples=n_samples)
    x_l1 = read_strain(RAW, "L1", gps, n_samples=n_samples)
    n0 = min(len(x_h1), len(x_l1))
    x_h1 = x_h1[:n0]
    x_l1 = x_l1[:n0]

    w_h1 = exact_phase_only_whiten(x_h1)
    w_l1 = exact_phase_only_whiten(x_l1)

    prof = lag_profile(w_h1, w_l1, fs=fs, max_ms=25.0)
    best_abs_lag_ms = float(prof["best_abs_lag_samples"] / fs * 1000.0)
    best_pos_lag_ms = float(prof["best_pos_lag_samples"] / fs * 1000.0)
    best_abs_corr = float(prof["best_abs_corr"])

    # Permutation by circular shifts > 1 sec to break short-lag alignment.
    min_shift = fs
    max_shift = n0 - fs
    null_max_abs = []
    for _ in range(n_perm):
        sh = int(rng.integers(min_shift, max_shift))
        y_shift = np.roll(w_l1, sh)
        p = lag_profile(w_h1, y_shift, fs=fs, max_ms=25.0)
        null_max_abs.append(abs(float(p["best_abs_corr"])))
    null_max_abs = np.array(null_max_abs, dtype=float)

    p_abs = float((np.sum(null_max_abs >= abs(best_abs_corr)) + 1.0) / (len(null_max_abs) + 1.0))
    near_light_travel = bool(abs(abs(best_abs_lag_ms) - 10.0) <= 2.0)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "gps": gps,
            "fs_hz": fs,
            "n_samples": n0,
            "n_perm": n_perm,
            "seed": seed,
            "whitening": "exact_phase_only_whiten",
            "lag_window_ms": [-25.0, 25.0],
            "null_model": "circular_shift_gt_1s",
        },
        "results": {
            "best_abs_lag_ms": best_abs_lag_ms,
            "best_abs_corr": best_abs_corr,
            "best_pos_lag_ms": best_pos_lag_ms,
            "best_pos_corr": float(prof["best_pos_corr"]),
            "corr_at_plus_10ms": float(
                prof["corr"][np.argmin(np.abs(prof["lags_samples"] - int(round(0.010 * fs))))]
            ),
            "near_light_travel_window_pm2ms": near_light_travel,
        },
        "null_significance": {
            "null_max_abs_corr_mean": float(np.mean(null_max_abs)),
            "null_max_abs_corr_std": float(np.std(null_max_abs)),
            "p_abs": p_abs,
        },
        "verdict": (
            "MICRO_TIMESHIFT_SIGNAL_SUPPORTED"
            if (p_abs <= 0.01 and near_light_travel)
            else "MICRO_TIMESHIFT_SIGNAL_NOT_SUPPORTED"
        ),
        "lag_profile_window_ms": {
            f"{(int(l)/fs)*1000.0:.2f} ms": float(c)
            for l, c in zip(prof["lags_samples"], prof["corr"])
        },
    }

    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

