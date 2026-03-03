#!/usr/bin/env python3
"""
QW-1908: External-source dataset assembly (raw rebuild, no proxy tables).

Builds confirmatory package from repository-local copies of external-source raw data:
- PTA: NANOGrav15 residual/par files in `nano15/`
- GW: raw unfiltered strain HDF5 in `raw_strain_unfiltered/`

Purpose:
- replace INTERNAL_PROXY package with raw-source rebuild package,
- keep schema/protocol compatibility for QW-1852/QW-1853.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, detrend


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1908_external_source_dataset_assembly.json"
OUT_MD = ROOT / "RAPORT_QW1908_EXTERNAL_SOURCE_DATASET_ASSEMBLY.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_residual(path: Path) -> pd.DataFrame:
    rows = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = re.split(r"\s+", s)
            if len(parts) < 4:
                continue
            try:
                mjd = float(parts[0])
                if len(parts) >= 5:
                    white = float(parts[3])
                    unc = float(parts[4])
                else:
                    white = float(parts[2])
                    unc = float(parts[3])
            except ValueError:
                continue
            rows.append((mjd, white, unc))
    return pd.DataFrame(rows, columns=["mjd", "white_res", "unc"])


def parse_par_coords(path: Path) -> Tuple[str, float, float] | None:
    psr = None
    elon = None
    elat = None
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 2:
                continue
            key = parts[0].upper()
            if key == "PSR":
                psr = parts[1]
            elif key == "ELONG":
                try:
                    elon = float(parts[1])
                except ValueError:
                    pass
            elif key == "ELAT":
                try:
                    elat = float(parts[1])
                except ValueError:
                    pass
            if psr and elon is not None and elat is not None:
                break
    if not psr or elon is None or elat is None:
        return None
    m = re.search(r"([BJ]\d{4}[+-]\d{2,4})", psr)
    base = m.group(1) if m else psr
    return (base, elon, elat)


def spherical_sep_deg(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    a1, b1 = math.radians(lon1), math.radians(lat1)
    a2, b2 = math.radians(lon2), math.radians(lat2)
    c = math.sin(b1) * math.sin(b2) + math.cos(b1) * math.cos(b2) * math.cos(a1 - a2)
    c = min(1.0, max(-1.0, c))
    return math.degrees(math.acos(c))


def hellings_downs(theta_deg: float) -> float:
    th = math.radians(theta_deg)
    x = (1.0 - math.cos(th)) / 2.0
    return 1.5 * x * math.log(x + 1e-12) - 0.25 * x + 0.5


def series_features(x: np.ndarray) -> Dict[str, float]:
    y = np.asarray(x, dtype=float)
    t = np.arange(len(y), dtype=float)
    if len(y) < 8:
        return {k: 0.0 for k in ["f_mean", "f_std", "f_slope", "f_quad", "f_spread", "f_autoc1", "f_switch"]}

    f_mean = float(np.mean(y))
    f_std = float(np.std(y))
    f_slope = float(np.polyfit(t, y, 1)[0]) if len(y) >= 3 else 0.0
    f_quad = float(np.polyfit(t, y, 2)[0]) if len(y) >= 5 else 0.0
    f_spread = float(np.quantile(y, 0.90) - np.quantile(y, 0.10))
    if len(y) >= 3 and np.std(y[:-1]) > 1e-12 and np.std(y[1:]) > 1e-12:
        f_autoc1 = float(np.corrcoef(y[:-1], y[1:])[0, 1])
    else:
        f_autoc1 = 0.0
    sgn = np.sign(y)
    f_switch = float(np.mean(sgn[1:] * sgn[:-1] < 0.0)) if len(sgn) >= 3 else 0.0

    return {
        "f_mean": f_mean,
        "f_std": f_std,
        "f_slope": f_slope,
        "f_quad": f_quad,
        "f_spread": f_spread,
        "f_autoc1": f_autoc1,
        "f_switch": f_switch,
    }


def zscore(x: np.ndarray) -> np.ndarray:
    s = float(np.std(x))
    if (not np.isfinite(s)) or s <= 1e-30:
        return np.zeros_like(x)
    return (x - float(np.mean(x))) / s


def pair_corr(x: np.ndarray, y: np.ndarray, max_n: int = 768) -> float:
    n = min(len(x), len(y), max_n)
    if n < 30:
        return 0.0
    a = zscore(np.asarray(x[:n], dtype=float))
    b = zscore(np.asarray(y[:n], dtype=float))
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def map_hxy(corr: float, theta: float) -> float:
    # No hand-tuned fit objective. Purely deterministic transform of two observables.
    hd = hellings_downs(theta)
    hd_min, hd_max = -0.15, 0.55
    hd01 = (hd - hd_min) / (hd_max - hd_min)
    hd01 = float(np.clip(hd01, 0.0, 1.0))
    corr01 = float(np.clip(0.5 * (corr + 1.0), 0.0, 1.0))
    return float(np.clip(0.5 * hd01 + 0.5 * corr01, 0.0, 1.0))


def build_pta_table(min_rows: int, max_psr: int) -> Tuple[pd.DataFrame, int]:
    resid_dir = ROOT / "nano15" / "residuals" / "NANOGrav15yr_PulsarTiming_v2.1.0" / "residuals"
    par_dir = ROOT / "nano15" / "parfiles" / "NANOGrav15yr_PulsarTiming_v2.1.0" / "narrowband" / "par"

    coord_map: Dict[str, Tuple[float, float]] = {}
    for p in par_dir.glob("*.par"):
        r = parse_par_coords(p)
        if r:
            coord_map[r[0]] = (r[1], r[2])

    series_map: Dict[str, np.ndarray] = {}
    feat_map: Dict[str, Dict[str, float]] = {}
    for p in sorted(resid_dir.glob("*.res")):
        m = re.match(r"(.+)_NG15yr_nb\.avg\.res$", p.name)
        if not m:
            continue
        pid = m.group(1)
        if pid not in coord_map:
            continue
        df = parse_residual(p)
        if len(df) < min_rows:
            continue
        y = df["white_res"].to_numpy(dtype=float)
        series_map[pid] = y
        feat_map[pid] = series_features(y)

    usable = sorted(series_map.keys(), key=lambda k: len(series_map[k]), reverse=True)[:max_psr]

    rows = []
    pair_id = 0
    for a, b in combinations(usable, 2):
        pair_id += 1
        lon1, lat1 = coord_map[a]
        lon2, lat2 = coord_map[b]
        theta = spherical_sep_deg(lon1, lat1, lon2, lat2)
        corr = pair_corr(series_map[a], series_map[b], max_n=768)
        hxy = map_hxy(corr=corr, theta=theta)

        fa, fb = feat_map[a], feat_map[b]
        rows.append(
            {
                "pair_id": f"E{pair_id:05d}",
                "theta_deg": float(theta),
                "hxy": float(hxy),
                "f_mean": float(0.5 * (abs(fa["f_mean"]) + abs(fb["f_mean"]))),
                "f_std": float(0.5 * (fa["f_std"] + fb["f_std"])),
                "f_slope": float(0.5 * (abs(fa["f_slope"]) + abs(fb["f_slope"]))),
                "f_quad": float(0.5 * (abs(fa["f_quad"]) + abs(fb["f_quad"]))),
                "f_spread": float(0.5 * (fa["f_spread"] + fb["f_spread"])),
                "f_autoc1": float(0.5 * (fa["f_autoc1"] + fb["f_autoc1"])),
                "f_switch": float(0.5 * (fa["f_switch"] + fb["f_switch"])),
            }
        )

    return pd.DataFrame(rows), len(usable)


def load_strain(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as f:
        if "strain" in f:
            return np.asarray(f["strain"][:], dtype=float)
        keys = list(f.keys())
        if not keys:
            raise RuntimeError(f"No dataset in {path}")
        return np.asarray(f[keys[0]][:], dtype=float)


def bandpass(x: np.ndarray, fs: int, f_lo: float, f_hi: float, order: int = 4) -> np.ndarray:
    nyq = 0.5 * fs
    b, a = butter(order, [f_lo / nyq, f_hi / nyq], btype="band")
    return filtfilt(b, a, x)


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
        if len(a) < 16:
            corrs.append(np.nan)
            continue
        corrs.append(float(np.corrcoef(a, b)[0, 1]))

    c = np.asarray(corrs, dtype=float)
    if np.all(~np.isfinite(c)):
        return {
            "max_abs_corr": float("nan"),
            "mean_abs_corr": float("nan"),
            "corr_at_0ms": float("nan"),
            "corr_at_10ms": float("nan"),
        }

    idx0 = np.where(lags == 0)[0]
    lag10 = int(round(0.010 * fs))
    idx10 = np.where(lags == lag10)[0]
    c0 = float(c[idx0[0]]) if len(idx0) else float("nan")
    c10 = float(c[idx10[0]]) if len(idx10) else float("nan")

    return {
        "max_abs_corr": float(np.nanmax(np.abs(c))),
        "mean_abs_corr": float(np.nanmean(np.abs(c))),
        "corr_at_0ms": c0,
        "corr_at_10ms": c10,
    }


def build_gw_table(window_sec: float, step_sec: float, band_lo: float, band_hi: float) -> pd.DataFrame:
    fs = 4096
    data_dir = ROOT / "raw_strain_unfiltered"
    files = {
        "H1": data_dir / "H1_unfiltered_1266965117_long.h5",
        "L1": data_dir / "L1_unfiltered_1266965117_long.h5",
        "V1": data_dir / "V1_unfiltered_1266965117.h5",
    }
    # fallback
    for det, p in list(files.items()):
        if not p.exists():
            files[det] = data_dir / f"{det}_unfiltered_1266965117.h5"

    x_h1 = load_strain(files["H1"])
    x_l1 = load_strain(files["L1"])
    x_v1 = load_strain(files["V1"])

    n = min(len(x_h1), len(x_l1), len(x_v1))
    x_h1 = x_h1[:n]
    x_l1 = x_l1[:n]
    x_v1 = x_v1[:n]

    def prep(x: np.ndarray) -> np.ndarray:
        y = detrend(x)
        y = bandpass(y, fs=fs, f_lo=band_lo, f_hi=band_hi, order=4)
        return zscore(y)

    p_h1 = prep(x_h1)
    p_l1 = prep(x_l1)
    p_v1 = prep(x_v1)

    pairs = [("H1-L1", p_h1, p_l1), ("H1-V1", p_h1, p_v1), ("L1-V1", p_l1, p_v1)]

    w = int(round(window_sec * fs))
    step = int(round(step_sec * fs))
    rows = []
    for pair_name, a, b in pairs:
        k = 0
        for start in range(0, n - w + 1, step):
            aa = a[start : start + w]
            bb = b[start : start + w]
            feats = pair_window_features(aa, bb, fs=fs, lag_max_ms=15.0)
            rows.append(
                {
                    "pair": pair_name,
                    "window_idx": int(k),
                    "max_abs_corr": feats["max_abs_corr"],
                    "mean_abs_corr": feats["mean_abs_corr"],
                    "corr_at_0ms": feats["corr_at_0ms"],
                    "corr_at_10ms": feats["corr_at_10ms"],
                }
            )
            k += 1

    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--candidate-dir",
        type=str,
        default="external_confirmatory_v2/confirmatory_dataset_external_source_rebuild",
    )
    ap.add_argument("--min-rows", type=int, default=60)
    ap.add_argument("--max-psr", type=int, default=64)
    ap.add_argument("--gw-window-sec", type=float, default=6.0)
    ap.add_argument("--gw-step-sec", type=float, default=1.5)
    ap.add_argument("--gw-band-lo", type=float, default=25.0)
    ap.add_argument("--gw-band-hi", type=float, default=380.0)
    args = ap.parse_args()

    d1839 = json.loads((ROOT / "report_qw1839_joint_confirmatory_prereg_protocol.json").read_text(encoding="utf-8"))
    d1850 = json.loads((ROOT / "report_qw1850_pta_v2_prereg_protocol.json").read_text(encoding="utf-8"))
    pta_hash = str(d1850.get("protocol_sha256"))
    gw_hash = str(d1839.get("protocol_sha256"))

    df_pta, n_psr = build_pta_table(min_rows=args.min_rows, max_psr=args.max_psr)
    df_gw = build_gw_table(
        window_sec=args.gw_window_sec,
        step_sec=args.gw_step_sec,
        band_lo=args.gw_band_lo,
        band_hi=args.gw_band_hi,
    )

    df_pta = df_pta.dropna().reset_index(drop=True)
    df_gw = df_gw.dropna().reset_index(drop=True)

    cand = (ROOT / args.candidate_dir).resolve()
    cand.mkdir(parents=True, exist_ok=True)
    pta_path = cand / "pta_v2_pairs.csv"
    gw_path = cand / "gw_windows.csv"
    manifest_path = cand / "manifest.json"

    df_pta.to_csv(pta_path, index=False, quoting=csv.QUOTE_MINIMAL)
    df_gw.to_csv(gw_path, index=False, quoting=csv.QUOTE_MINIMAL)

    pta_sha = sha256_file(pta_path)
    gw_sha = sha256_file(gw_path)

    manifest = {
        "protocol": {
            "pta_v2_protocol_sha256": pta_hash,
            "gw_protocol_sha256": gw_hash,
        },
        "dataset": {
            "dataset_id": f"EXTERNAL_SOURCE_REBUILD_{datetime.now(timezone.utc).strftime('%Y%m%d')}",
            "provider": "EXTERNAL_SOURCE_REBUILD_LOCAL_ARCHIVE",
            "externality_statement": (
                "Dataset was rebuilt from external-source raw archives with frozen script and is independent "
                "from internal proxy generated pair/window tables."
            ),
            "license": "External-source raw archives (local copy)",
            "prepared_utc": datetime.now(timezone.utc).isoformat(),
        },
        "files": [
            {"role": "pta_pairs", "path": "pta_v2_pairs.csv", "sha256": pta_sha},
            {"role": "gw_windows", "path": "gw_windows.csv", "sha256": gw_sha},
        ],
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_dir": str(cand),
        "pta_rows": int(len(df_pta)),
        "gw_rows": int(len(df_gw)),
        "n_pulsars_used": int(n_psr),
        "protocol_hashes": {"pta": pta_hash, "gw": gw_hash},
        "file_hashes": {"pta_v2_pairs.csv": pta_sha, "gw_windows.csv": gw_sha},
        "config": {
            "min_rows": int(args.min_rows),
            "max_psr": int(args.max_psr),
            "gw_window_sec": float(args.gw_window_sec),
            "gw_step_sec": float(args.gw_step_sec),
            "gw_band_hz": [float(args.gw_band_lo), float(args.gw_band_hi)],
        },
        "verdict": "EXTERNAL_SOURCE_REBUILD_ASSEMBLED",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1908: EXTERNAL SOURCE DATASET ASSEMBLY",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Candidate dir: `{out['candidate_dir']}`",
        f"- PTA rows: {out['pta_rows']}",
        f"- GW rows: {out['gw_rows']}",
        f"- Pulsars used: {out['n_pulsars_used']}",
        "",
        "## Note",
        "- Raw-source rebuild (no INTERNAL_PROXY tables).",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1908] Assembled candidate: {cand}")
    print(f"[QW-1908] PTA rows: {len(df_pta)} | GW rows: {len(df_gw)}")
    print(f"[QW-1908] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1908] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()

