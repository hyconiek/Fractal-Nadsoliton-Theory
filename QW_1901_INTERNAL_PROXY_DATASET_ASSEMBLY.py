#!/usr/bin/env python3
"""
QW-1901: Internal proxy dataset assembly for frozen-pipeline rehearsal.

Builds a candidate package compatible with QW-1852/1853 schema using local files:
- PTA pairs from nano15 residual/par files
- GW windows from gw1831_window_features.csv

Important: marks dataset as INTERNAL PROXY (not independent confirmatory evidence).
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
CAND_DIR = ROOT / "external_confirmatory_v2" / "confirmatory_dataset_internal_proxy"
OUT_JSON = ROOT / "report_qw1901_internal_proxy_dataset_assembly.json"
OUT_MD = ROOT / "RAPORT_QW1901_INTERNAL_PROXY_DATASET_ASSEMBLY.md"


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
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 4:
                continue
            try:
                mjd = float(parts[0])
                # Two formats exist in local residual files:
                # 5 cols: mjd, freq, residual, white_residual, uncertainty
                # 4 cols: mjd, freq, residual, uncertainty
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
    return (psr, elon, elat)


def base_psr_id(psr: str) -> str:
    m = re.search(r"([BJ]\d{4}[+-]\d{2,4})", psr)
    return m.group(1) if m else psr


def spherical_sep_deg(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    a1, b1 = math.radians(lon1), math.radians(lat1)
    a2, b2 = math.radians(lon2), math.radians(lat2)
    c = math.sin(b1) * math.sin(b2) + math.cos(b1) * math.cos(b2) * math.cos(a1 - a2)
    c = min(1.0, max(-1.0, c))
    return math.degrees(math.acos(c))


def series_features(x: np.ndarray) -> Dict[str, float]:
    if len(x) < 6:
        return {
            "f_mean": 0.0,
            "f_std": 0.0,
            "f_slope": 0.0,
            "f_quad": 0.0,
            "f_spread": 0.0,
            "f_autoc1": 0.0,
            "f_switch": 0.0,
        }

    y = np.asarray(x, dtype=float)
    t = np.arange(len(y), dtype=float)

    f_mean = float(np.mean(y))
    f_std = float(np.std(y))

    if len(y) >= 3:
        p1 = np.polyfit(t, y, 1)
        f_slope = float(p1[0])
    else:
        f_slope = 0.0

    if len(y) >= 5:
        p2 = np.polyfit(t, y, 2)
        f_quad = float(p2[0])
    else:
        f_quad = 0.0

    f_spread = float(np.quantile(y, 0.90) - np.quantile(y, 0.10))

    if len(y) >= 3 and np.std(y[:-1]) > 1e-12 and np.std(y[1:]) > 1e-12:
        f_autoc1 = float(np.corrcoef(y[:-1], y[1:])[0, 1])
    else:
        f_autoc1 = 0.0

    sgn = np.sign(y)
    if len(sgn) >= 3:
        f_switch = float(np.mean(sgn[1:] * sgn[:-1] < 0.0))
    else:
        f_switch = 0.0

    return {
        "f_mean": f_mean,
        "f_std": f_std,
        "f_slope": f_slope,
        "f_quad": f_quad,
        "f_spread": f_spread,
        "f_autoc1": f_autoc1,
        "f_switch": f_switch,
    }


def pair_corr(x: np.ndarray, y: np.ndarray) -> float:
    n = min(len(x), len(y), 256)
    if n < 12:
        return 0.0
    a = np.asarray(x[:n], dtype=float)
    b = np.asarray(y[:n], dtype=float)
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def main() -> None:
    resid_dir = ROOT / "nano15" / "residuals" / "NANOGrav15yr_PulsarTiming_v2.1.0" / "residuals"
    par_dir = ROOT / "nano15" / "parfiles" / "NANOGrav15yr_PulsarTiming_v2.1.0" / "narrowband" / "par"
    gw_src = ROOT / "gw1831_window_features.csv"

    d1839 = json.loads((ROOT / "report_qw1839_joint_confirmatory_prereg_protocol.json").read_text(encoding="utf-8"))
    d1850 = json.loads((ROOT / "report_qw1850_pta_v2_prereg_protocol.json").read_text(encoding="utf-8"))

    pta_hash = str(d1850.get("protocol_sha256"))
    gw_hash = str(d1839.get("protocol_sha256"))

    # Coordinates map from par files.
    coord_map: Dict[str, Tuple[float, float]] = {}
    for p in par_dir.glob("*.par"):
        r = parse_par_coords(p)
        if not r:
            continue
        psr, elon, elat = r
        coord_map[base_psr_id(psr)] = (elon, elat)

    # Residual series map.
    series_map: Dict[str, np.ndarray] = {}
    feat_map: Dict[str, Dict[str, float]] = {}

    for p in resid_dir.glob("*.res"):
        m = re.match(r"(.+)_NG15yr_nb\.avg\.res$", p.name)
        if not m:
            continue
        pid = m.group(1)
        df = parse_residual(p)
        if len(df) < 24:
            continue
        y = df["white_res"].to_numpy(dtype=float)
        series_map[pid] = y
        feat_map[pid] = series_features(y)

    usable = [pid for pid in series_map.keys() if pid in coord_map]
    usable = sorted(usable, key=lambda x: len(series_map[x]), reverse=True)

    # Keep top 24 pulsars -> 276 pairs (well above min 80).
    usable = usable[:24]

    pta_rows = []
    pair_id = 0
    for a, b in combinations(usable, 2):
        pair_id += 1
        lon1, lat1 = coord_map[a]
        lon2, lat2 = coord_map[b]
        theta = spherical_sep_deg(lon1, lat1, lon2, lat2)

        corr = pair_corr(series_map[a], series_map[b])
        # proxy Hxy in plausible range [0,1]
        hxy = float(np.clip(0.52 + 0.22 * corr + 0.18 * (1.0 - theta / 180.0), 0.0, 1.0))

        fa, fb = feat_map[a], feat_map[b]
        row = {
            "pair_id": f"P{pair_id:04d}",
            "theta_deg": theta,
            "hxy": hxy,
            "f_mean": 0.5 * (abs(fa["f_mean"]) + abs(fb["f_mean"])),
            "f_std": 0.5 * (fa["f_std"] + fb["f_std"]),
            "f_slope": 0.5 * (abs(fa["f_slope"]) + abs(fb["f_slope"])),
            "f_quad": 0.5 * (abs(fa["f_quad"]) + abs(fb["f_quad"])),
            "f_spread": 0.5 * (fa["f_spread"] + fb["f_spread"]),
            "f_autoc1": 0.5 * (fa["f_autoc1"] + fb["f_autoc1"]),
            "f_switch": 0.5 * (fa["f_switch"] + fb["f_switch"]),
        }
        pta_rows.append(row)

    df_pta = pd.DataFrame(pta_rows)

    # GW windows from existing event-window features (schema subset).
    df_gw_src = pd.read_csv(gw_src)
    df_gw = df_gw_src[["pair", "window_idx", "max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]].copy()

    CAND_DIR.mkdir(parents=True, exist_ok=True)
    pta_path = CAND_DIR / "pta_v2_pairs.csv"
    gw_path = CAND_DIR / "gw_windows.csv"
    manifest_path = CAND_DIR / "manifest.json"

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
            "dataset_id": "INTERNAL_PROXY_NANO15_GW1831_20260303",
            "provider": "INTERNAL_PROXY",
            "externality_statement": (
                "Dataset assembled from repository-local files and is NOT independent from prior in-repo analyses; "
                "it is for frozen-pipeline rehearsal only, not external confirmatory evidence."
            ),
            "license": "Repository-local mixed",
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
        "candidate_dir": str(CAND_DIR),
        "pta_rows": int(len(df_pta)),
        "gw_rows": int(len(df_gw)),
        "n_pulsars_used": int(len(usable)),
        "protocol_hashes": {
            "pta": pta_hash,
            "gw": gw_hash,
        },
        "file_hashes": {
            "pta_v2_pairs.csv": pta_sha,
            "gw_windows.csv": gw_sha,
        },
        "verdict": "INTERNAL_PROXY_DATASET_ASSEMBLED",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1901: INTERNAL PROXY DATASET ASSEMBLY",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Candidate dir: `{out['candidate_dir']}`",
        "",
        "## Rows",
        f"- PTA pairs: {out['pta_rows']}",
        f"- GW windows: {out['gw_rows']}",
        f"- Pulsars used: {out['n_pulsars_used']}",
        "",
        "## Note",
        "- This package is INTERNAL PROXY only (not independent external confirmatory dataset).",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1901] Assembled candidate: {CAND_DIR}")
    print(f"[QW-1901] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1901] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
