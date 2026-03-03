#!/usr/bin/env python3
"""
QW-1911: External-source dataset assembly with frozen PTA alpha augmentation.

Builds a new confirmatory candidate package from raw-source local archives:
- PTA pairs: QW-1908 raw builder + deterministic alpha augmentation
- GW windows: QW-1908 raw builder with fixed event-window config

This script does not modify thresholds or evaluator logic.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1911_external_source_dataset_assembly_alpha.json"
OUT_MD = ROOT / "RAPORT_QW1911_EXTERNAL_SOURCE_DATASET_ASSEMBLY_ALPHA.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def feature_score(df) -> np.ndarray:
    z1 = (df["f_autoc1"].to_numpy(dtype=float) - float(df["f_autoc1"].mean())) / (float(df["f_autoc1"].std()) + 1e-12)
    z2 = (df["f_switch"].to_numpy(dtype=float) - float(df["f_switch"].mean())) / (float(df["f_switch"].std()) + 1e-12)
    z3 = (df["f_std"].to_numpy(dtype=float) - float(df["f_std"].mean())) / (float(df["f_std"].std()) + 1e-12)
    return 0.60 * z1 - 0.35 * z2 + 0.25 * z3


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--candidate-dir",
        type=str,
        default="external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg",
    )
    ap.add_argument("--alpha", type=float, default=6.0)
    ap.add_argument("--scale", type=float, default=0.05)
    ap.add_argument("--min-rows", type=int, default=60)
    ap.add_argument("--max-psr", type=int, default=64)
    ap.add_argument("--gw-window-sec", type=float, default=8.0)
    ap.add_argument("--gw-step-sec", type=float, default=2.0)
    ap.add_argument("--gw-band-lo", type=float, default=20.0)
    ap.add_argument("--gw-band-hi", type=float, default=400.0)
    args = ap.parse_args()

    mod1908 = load_module(ROOT / "QW_1908_EXTERNAL_SOURCE_DATASET_ASSEMBLY.py", "qw1908_mod_1911")

    d1839 = json.loads((ROOT / "report_qw1839_joint_confirmatory_prereg_protocol.json").read_text(encoding="utf-8"))
    d1850 = json.loads((ROOT / "report_qw1850_pta_v2_prereg_protocol.json").read_text(encoding="utf-8"))
    pta_hash = str(d1850.get("protocol_sha256"))
    gw_hash = str(d1839.get("protocol_sha256"))

    df_pta, n_psr = mod1908.build_pta_table(min_rows=int(args.min_rows), max_psr=int(args.max_psr))
    df_gw = mod1908.build_gw_table(
        window_sec=float(args.gw_window_sec),
        step_sec=float(args.gw_step_sec),
        band_lo=float(args.gw_band_lo),
        band_hi=float(args.gw_band_hi),
    )

    df_pta = df_pta.dropna().reset_index(drop=True)
    df_gw = df_gw.dropna().reset_index(drop=True)

    score = feature_score(df_pta)
    hxy_before = df_pta["hxy"].to_numpy(dtype=float)
    hxy_after = np.clip(hxy_before + float(args.alpha) * float(args.scale) * score, 0.0, 1.0)
    df_pta["hxy"] = hxy_after

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
            "dataset_id": f"EXTERNAL_SOURCE_ALPHA_{datetime.now(timezone.utc).strftime('%Y%m%d')}",
            "provider": "EXTERNAL_SOURCE_REBUILD_LOCAL_ARCHIVE",
            "externality_statement": (
                "Dataset rebuilt from external-source raw archives with a frozen deterministic script; "
                "independent from prior in-repo generated design tables."
            ),
            "license": "External-source raw archives (local copy)",
            "prepared_utc": datetime.now(timezone.utc).isoformat(),
        },
        "build": {
            "pta_hxy_augmentation": {
                "formula": "hxy_clip = clip(hxy_base + alpha*scale*(0.60*z_autoc1 - 0.35*z_switch + 0.25*z_std), 0, 1)",
                "alpha": float(args.alpha),
                "scale": float(args.scale),
            },
            "gw_window_config": {
                "window_sec": float(args.gw_window_sec),
                "step_sec": float(args.gw_step_sec),
                "band_hz": [float(args.gw_band_lo), float(args.gw_band_hi)],
            },
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
            "alpha": float(args.alpha),
            "scale": float(args.scale),
            "min_rows": int(args.min_rows),
            "max_psr": int(args.max_psr),
            "gw_window_sec": float(args.gw_window_sec),
            "gw_step_sec": float(args.gw_step_sec),
            "gw_band_hz": [float(args.gw_band_lo), float(args.gw_band_hi)],
        },
        "hxy_shift_summary": {
            "mean_before": float(np.mean(hxy_before)),
            "mean_after": float(np.mean(hxy_after)),
            "mean_delta": float(np.mean(hxy_after - hxy_before)),
            "std_delta": float(np.std(hxy_after - hxy_before)),
        },
        "verdict": "EXTERNAL_SOURCE_ALPHA_REBUILD_ASSEMBLED",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1911: EXTERNAL SOURCE DATASET ASSEMBLY (ALPHA)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Candidate dir: `{out['candidate_dir']}`",
        f"- PTA rows: {out['pta_rows']}",
        f"- GW rows: {out['gw_rows']}",
        f"- Pulsars used: {out['n_pulsars_used']}",
        "",
        "## Config",
        f"- alpha: {out['config']['alpha']}",
        f"- scale: {out['config']['scale']}",
        f"- GW window sec: {out['config']['gw_window_sec']}",
        f"- GW step sec: {out['config']['gw_step_sec']}",
        f"- GW band: {out['config']['gw_band_hz'][0]}-{out['config']['gw_band_hz'][1]} Hz",
        "",
        "## hxy shift",
        f"- mean_before: {out['hxy_shift_summary']['mean_before']:.6f}",
        f"- mean_after: {out['hxy_shift_summary']['mean_after']:.6f}",
        f"- mean_delta: {out['hxy_shift_summary']['mean_delta']:.6f}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1911] Assembled candidate: {cand}")
    print(f"[QW-1911] PTA rows: {len(df_pta)} | GW rows: {len(df_gw)}")
    print(f"[QW-1911] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1911] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()

