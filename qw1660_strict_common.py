#!/usr/bin/env python3
"""
Common strict utilities for QW-1660 GW raw-strain methodology repairs.

Design goals:
- One canonical estimator implementation (Cross-MF-DFA q=0)
- One canonical reader for cached raw strain
- Deterministic, auditable preprocessing helpers
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
from scipy.signal import detrend
from gwpy.timeseries import TimeSeries


def _try_read_h5(path: Path) -> np.ndarray | None:
    try:
        with h5py.File(path, "r") as f:
            if "strain" in f:
                return np.asarray(f["strain"][:], dtype=float)
            key = list(f.keys())[0]
            return np.asarray(f[key][:], dtype=float)
    except Exception:
        return None


def _fetch_and_cache_raw(
    raw_dir: Path,
    det: str,
    gps: int,
    n_samples: int | None,
    fs: int,
) -> np.ndarray:
    need = n_samples if n_samples is not None else int(512 * fs)
    duration = max(512, int(np.ceil(need / fs)) + 8)
    try:
        ts = TimeSeries.fetch_open_data(det, gps, gps + duration, verbose=False, cache=True)
    except Exception as exc:
        raise RuntimeError(
            "Failed to fetch GWOSC data. Cached raw files are missing/invalid "
            "(often Git LFS pointers) and network access to gwosc.org is unavailable. "
            "Provide valid local HDF5 strain files in raw_strain_unfiltered/ or run in "
            "an environment with GWOSC connectivity."
        ) from exc
    if ts.sample_rate.value > fs:
        ts = ts.resample(fs)
    ts = ts.notch(60).notch(120).notch(180)
    x = np.asarray(ts.value, dtype=float)
    if n_samples is not None:
        x = x[:n_samples]
    out_path = raw_dir / f"{det}_unfiltered_{gps}.h5"
    with h5py.File(out_path, "w") as f:
        f.create_dataset("strain", data=x)
    return x


def read_strain(
    raw_dir: Path,
    det: str,
    gps: int,
    n_samples: int | None = None,
    fs: int = 4096,
    fetch_if_invalid: bool = True,
) -> np.ndarray:
    """Read cached raw strain with fallback naming conventions."""
    p1 = raw_dir / f"{det}_unfiltered_{gps}.h5"
    p2 = raw_dir / f"{det}_unfiltered.h5"
    path = p1 if p1.exists() else p2
    if not path.exists():
        if not fetch_if_invalid:
            raise FileNotFoundError(f"Missing cached strain file: {p1} (or fallback {p2})")
        x = _fetch_and_cache_raw(raw_dir, det, gps, n_samples=n_samples, fs=fs)
        return detrend(np.asarray(x, dtype=float))

    x = _try_read_h5(path)
    if x is None:
        if not fetch_if_invalid:
            raise OSError(f"Invalid/non-HDF5 cached file: {path}")
        # Typical case when Git LFS pointer placeholders are present.
        x = _fetch_and_cache_raw(raw_dir, det, gps, n_samples=n_samples, fs=fs)

    if n_samples is not None:
        x = x[:n_samples]
    return detrend(np.asarray(x, dtype=float))


def scales_for_n(n: int, n_scales: int = 15, min_scale_pow10: int = 3) -> np.ndarray:
    """Canonical scale grid for MF-DFA."""
    if n < 10 ** min_scale_pow10 * 4:
        raise ValueError(f"n={n} too small for min scale 10^{min_scale_pow10}.")
    return np.unique(np.logspace(min_scale_pow10, np.log10(n // 4), n_scales).astype(int))


def cross_mfdfa_q0(x: np.ndarray, y: np.ndarray, scales: np.ndarray) -> float:
    """
    Canonical Cross-MF-DFA q=0 estimator used across strict runs.
    """
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


def phase_randomized_signal(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """PSD-preserving phase randomization."""
    xf = np.fft.rfft(x)
    amp = np.abs(xf)
    ph = rng.uniform(0.0, 2.0 * np.pi, size=len(xf))
    ph[0] = 0.0
    if len(x) % 2 == 0:
        ph[-1] = 0.0
    out = np.fft.irfft(amp * np.exp(1j * ph), n=len(x))
    return detrend(out)


def exact_phase_only_whiten(x: np.ndarray) -> np.ndarray:
    """
    Phase-only whitening:
    - removes amplitude envelope by setting |X(f)|=1
    - preserves Fourier phase
    """
    xf = np.fft.rfft(x)
    ph = np.angle(xf)
    xw = np.exp(1j * ph)
    xw[0] = 0.0
    if len(x) % 2 == 0:
        xw[-1] = np.real(xw[-1])
    out = np.fft.irfft(xw, n=len(x))
    return detrend(out)


def safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    """Numerically stable Pearson correlation for large vectors."""
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


def empirical_p_lower(obs: float, null: np.ndarray) -> float:
    return float((np.sum(null <= obs) + 1.0) / (len(null) + 1.0))


def empirical_p_upper(obs: float, null: np.ndarray) -> float:
    return float((np.sum(null >= obs) + 1.0) / (len(null) + 1.0))
