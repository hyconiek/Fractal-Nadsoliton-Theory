#!/usr/bin/env python3
"""
QW-1725: Strict reanalysis of GW cross-Hurst anomaly.

Pre-registered protocol:
1) Fixed estimator: Cross-MF-DFA q=0 (same formula in all tests).
2) Fixed data source: cached raw strain files (no refetch).
3) Fixed controls:
   - Null A: independent phase randomization preserving PSD amplitudes.
   - Null B: circular time shift (>1s) breaking alignment.
4) Additional diagnostics:
   - sample-length sensitivity,
   - windowed stability,
   - +/-25ms lag correlation profile.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
from scipy.signal import detrend


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw_strain_unfiltered"
OUT_JSON = ROOT / "report_qw1725_gw_strict_cross_hurst_reanalysis.json"
OUT_MD = ROOT / "RAPORT_QW1725_GW_STRICT_CROSS_HURST_REANALYSIS.md"


def read_strain(path: Path, n: int | None = None) -> np.ndarray:
    with h5py.File(path, "r") as f:
        if "strain" in f:
            x = f["strain"][:]
        else:
            key = list(f.keys())[0]
            x = f[key][:]
    if n is not None:
        x = x[:n]
    return detrend(np.asarray(x, dtype=float))


def scales_for_n(n: int) -> np.ndarray:
    return np.unique(np.logspace(3, np.log10(n // 4), 15).astype(int))


def cross_mfdfa_q0(x: np.ndarray, y: np.ndarray, scales: np.ndarray) -> float:
    z = np.cumsum((x - np.mean(x)) * (y - np.mean(y)))
    n = len(z)
    vals = []
    used = []
    for s in scales:
        ns = n // s
        if ns == 0:
            continue
        t = np.arange(s)
        rms = []
        for i in range(ns):
            seg = z[i * s : (i + 1) * s]
            p = np.polyfit(t, seg, 1)
            trend = np.polyval(p, t)
            var = np.mean((seg - trend) ** 2)
            if var > 0:
                rms.append(var)
        if rms:
            r = np.array(rms, dtype=float)
            vals.append(np.exp(0.5 * np.mean(np.log(r + 1e-300))))
            used.append(s)
    if len(vals) < 3:
        return float("nan")
    return float(np.polyfit(np.log(used), np.log(vals), 1)[0])


def empirical_p_lower(obs: float, null: np.ndarray) -> float:
    return float((np.sum(null <= obs) + 1.0) / (len(null) + 1.0))


def empirical_p_upper(obs: float, null: np.ndarray) -> float:
    return float((np.sum(null >= obs) + 1.0) / (len(null) + 1.0))


def phase_randomized_signal(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    xf = np.fft.rfft(x)
    amp = np.abs(xf)
    ph = rng.uniform(0.0, 2.0 * np.pi, size=len(xf))
    ph[0] = 0.0
    if len(x) % 2 == 0:
        ph[-1] = 0.0
    out = np.fft.irfft(amp * np.exp(1j * ph), n=len(x))
    return detrend(out)


def lag_profile_corr(x: np.ndarray, y: np.ndarray, fs: int, max_ms: float = 25.0) -> dict:
    def safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
        a = np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0)
        b = np.nan_to_num(b, nan=0.0, posinf=0.0, neginf=0.0)
        a0 = a - np.mean(a)
        b0 = b - np.mean(b)
        da = float(np.sqrt(np.mean(a0 * a0)))
        db = float(np.sqrt(np.mean(b0 * b0)))
        if (not np.isfinite(da)) or (not np.isfinite(db)) or da < 1e-30 or db < 1e-30:
            return 0.0
        c = float(np.mean(a0 * b0) / (da * db))
        if not np.isfinite(c):
            return 0.0
        return c

    max_lag = int(round(max_ms * fs / 1000.0))
    lags = np.arange(-max_lag, max_lag + 1)
    center_start = max_lag
    center_end = len(y) - max_lag
    y0 = y[center_start:center_end]
    corr = []
    for lag in lags:
        xs = x[center_start + lag : center_end + lag]
        c = safe_pearson(xs, y0)
        corr.append(c)
    corr = np.array(corr, dtype=float)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    idx_abs = int(np.argmax(np.abs(corr)))
    idx_pos = int(np.argmax(corr))
    lag_ms_abs = float(lags[idx_abs] / fs * 1000.0)
    lag_ms_pos = float(lags[idx_pos] / fs * 1000.0)
    val_10ms = float(corr[np.argmin(np.abs(lags - int(round(0.010 * fs))))])
    return {
        "best_abs_lag_ms": lag_ms_abs,
        "best_abs_corr": float(corr[idx_abs]),
        "best_pos_lag_ms": lag_ms_pos,
        "best_pos_corr": float(corr[idx_pos]),
        "corr_at_plus_10ms": val_10ms,
        "window": {f"{lags[i]/fs*1000.0:.2f} ms": float(corr[i]) for i in range(len(lags))},
    }


def main() -> None:
    rng = np.random.default_rng(1725)
    fs = 4096
    gps = 1266965117

    p_h1 = RAW / f"H1_unfiltered_{gps}.h5"
    p_l1 = RAW / f"L1_unfiltered_{gps}.h5"
    p_v1 = RAW / f"V1_unfiltered_{gps}.h5"
    if not (p_h1.exists() and p_l1.exists()):
        raise FileNotFoundError("Brak danych H1/L1 w raw_strain_unfiltered.")

    x_h1_full = read_strain(p_h1)
    x_l1_full = read_strain(p_l1)
    n_full = min(len(x_h1_full), len(x_l1_full))
    x_h1_full = x_h1_full[:n_full]
    x_l1_full = x_l1_full[:n_full]

    # Length scan
    length_scan_rows = []
    length_list = [131072, 262144, 524288, 1048576, 2097152]
    for n in length_list:
        if n > n_full:
            continue
        h = cross_mfdfa_q0(x_h1_full[:n], x_l1_full[:n], scales_for_n(n))
        length_scan_rows.append({"n_samples": int(n), "h_cross": float(h)})

    # Baseline analysis at fixed N (pre-registered)
    n0 = 524288 if n_full >= 524288 else n_full
    x_h1 = x_h1_full[:n0]
    x_l1 = x_l1_full[:n0]
    scales0 = scales_for_n(n0)
    h_obs = cross_mfdfa_q0(x_h1, x_l1, scales0)

    # Windowed stability (non-overlapping)
    win_rows = []
    n_windows = n_full // n0
    for k in range(n_windows):
        a = k * n0
        b = (k + 1) * n0
        h_k = cross_mfdfa_q0(x_h1_full[a:b], x_l1_full[a:b], scales0)
        win_rows.append({"window": k, "start": a, "stop": b, "h_cross": float(h_k)})

    # Null A: PSD-preserving phase randomization
    n_null = 120
    null_phase = []
    for _ in range(n_null):
        s1 = phase_randomized_signal(x_h1, rng)
        s2 = phase_randomized_signal(x_l1, rng)
        null_phase.append(cross_mfdfa_q0(s1, s2, scales0))
    null_phase = np.array(null_phase, dtype=float)

    # Null B: circular time shifts (>1s)
    min_shift = fs
    max_shift = n0 - fs
    null_shift = []
    for _ in range(n_null):
        sh = int(rng.integers(min_shift, max_shift))
        y_shift = np.roll(x_l1, sh)
        null_shift.append(cross_mfdfa_q0(x_h1, y_shift, scales0))
    null_shift = np.array(null_shift, dtype=float)

    # Detector-pair cross-check (H1-V1) with same estimator
    hv = None
    if p_v1.exists():
        x_v1 = read_strain(p_v1, n=n0)
        hv = cross_mfdfa_q0(x_h1[: len(x_v1)], x_v1[: len(x_h1)], scales_for_n(len(x_v1)))

    # Lag profile
    lag = lag_profile_corr(x_h1, x_l1, fs=fs, max_ms=25.0)

    # Significance
    p_lower_phase = empirical_p_lower(h_obs, null_phase)
    p_lower_shift = empirical_p_lower(h_obs, null_shift)
    z_phase = float((h_obs - np.mean(null_phase)) / (np.std(null_phase) + 1e-12))
    z_shift = float((h_obs - np.mean(null_shift)) / (np.std(null_shift) + 1e-12))

    # Robustness summary
    h_vals = np.array([r["h_cross"] for r in length_scan_rows], dtype=float)
    length_spread = float(np.max(h_vals) - np.min(h_vals)) if len(h_vals) else float("nan")
    win_vals = np.array([r["h_cross"] for r in win_rows], dtype=float) if win_rows else np.array([])
    win_std = float(np.std(win_vals)) if len(win_vals) else float("nan")

    thresholds = {
        "null_p_max": 0.01,
        "length_spread_max": 0.08,
        "window_std_max": 0.05,
        "lag_pos_corr_at_10ms_min": 0.02,
    }
    robust_anomaly = (
        p_lower_phase <= thresholds["null_p_max"]
        and p_lower_shift <= thresholds["null_p_max"]
        and length_spread <= thresholds["length_spread_max"]
        and (np.isnan(win_std) or win_std <= thresholds["window_std_max"])
        and lag["corr_at_plus_10ms"] >= thresholds["lag_pos_corr_at_10ms_min"]
    )
    verdict = "GW_CROSS_HURST_ANOMALY_ROBUST" if robust_anomaly else "GW_CROSS_HURST_ANOMALY_NOT_ROBUST"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "gps": gps,
            "fs_hz": fs,
            "baseline_n_samples": n0,
            "n_null_trials_per_model": n_null,
            "seed": 1725,
        },
        "length_scan": length_scan_rows,
        "baseline": {
            "h_obs_h1_l1": float(h_obs),
            "h_h1_v1": None if hv is None else float(hv),
            "window_rows": win_rows,
        },
        "null_phase_randomized": {
            "mean": float(np.mean(null_phase)),
            "std": float(np.std(null_phase)),
            "p_lower": p_lower_phase,
            "z_score": z_phase,
        },
        "null_circular_shift": {
            "mean": float(np.mean(null_shift)),
            "std": float(np.std(null_shift)),
            "p_lower": p_lower_shift,
            "z_score": z_shift,
        },
        "lag_profile": lag,
        "stability": {
            "length_spread": length_spread,
            "window_std": win_std,
        },
        "thresholds": thresholds,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1725: GW STRICT CROSS-HURST REANALYSIS",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## Baseline",
        f"- H_cross(H1,L1) @ N={n0}: {h_obs:.4f}",
        f"- H_cross(H1,V1) @ N={n0}: {hv if hv is not None else 'n/a'}",
        "",
        "## Null models",
        f"- Phase-randomized: mean={np.mean(null_phase):.4f}, std={np.std(null_phase):.4f}, p_lower={p_lower_phase:.4g}, z={z_phase:.2f}",
        f"- Circular-shift: mean={np.mean(null_shift):.4f}, std={np.std(null_shift):.4f}, p_lower={p_lower_shift:.4g}, z={z_shift:.2f}",
        "",
        "## Stability",
        f"- length spread: {length_spread:.4f}",
        f"- window std: {win_std:.4f}",
        "",
        "## Lag test (+/-25 ms)",
        f"- corr(+10 ms): {lag['corr_at_plus_10ms']:.5f}",
        f"- best abs lag: {lag['best_abs_lag_ms']:.2f} ms (corr={lag['best_abs_corr']:.5f})",
        "",
        "## Artefakty",
        f"- JSON szczegolowy: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1725] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1725] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
