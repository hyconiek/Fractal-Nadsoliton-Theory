#!/usr/bin/env python3
"""
Ensure QW-1660 raw inputs are present and valid in raw_strain_unfiltered/.

This script does not require large binaries in git:
- if local HDF5 exists and is valid -> reuse
- if missing/invalid (e.g., LFS pointer text) -> fetch from GWOSC via gwpy
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import h5py

from qw1660_strict_common import read_strain


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw_strain_unfiltered"


def is_valid_h5(path: Path) -> bool:
    try:
        with h5py.File(path, "r"):
            return True
    except Exception:
        return False


def file_info(path: Path) -> dict:
    return {
        "exists": path.exists(),
        "size_bytes": (path.stat().st_size if path.exists() else 0),
        "valid_hdf5": (is_valid_h5(path) if path.exists() else False),
    }


def main() -> None:
    gps = int(os.getenv("QW1660_GPS", "1266965117"))
    n_samples = int(os.getenv("QW1660_STRICT_N_SAMPLES", "131072"))
    fs = int(os.getenv("QW1660_FS_HZ", "4096"))
    dets = ["H1", "L1"]

    RAW.mkdir(parents=True, exist_ok=True)

    before = {d: file_info(RAW / f"{d}_unfiltered_{gps}.h5") for d in dets}

    loaded = {}
    for d in dets:
        x = read_strain(RAW, d, gps, n_samples=n_samples, fs=fs, fetch_if_invalid=True)
        loaded[d] = {
            "n_samples_loaded": int(len(x)),
            "mean": float(x.mean()),
            "std": float(x.std()),
        }

    after = {d: file_info(RAW / f"{d}_unfiltered_{gps}.h5") for d in dets}

    out = {
        "gps": gps,
        "n_samples_requested": n_samples,
        "fs_hz": fs,
        "before": before,
        "after": after,
        "loaded": loaded,
        "status": "QW1660_RAW_INPUTS_READY",
    }
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

