#!/usr/bin/env python3
"""
QW-1930: True external beta-channel autocollector.

Builds external_confirmatory_v2/beta_channel_true_external package directly from:
1) NANOGrav 15yr raw timing archive (local tar.gz),
2) GWOSC GWTC event API (online with local cache fallback).

No internal proxy/rebuild pair tables are used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple
from urllib.request import urlopen

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1930_true_external_beta_channel_autocollector.json"
OUT_MD = ROOT / "RAPORT_QW1930_TRUE_EXTERNAL_BETA_CHANNEL_AUTOCOLLECTOR.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_u64(text: str) -> int:
    d = hashlib.sha256(text.encode("utf-8")).digest()
    return int.from_bytes(d[:8], "big", signed=False)


def parse_flag_map(tokens: List[str], start: int) -> Dict[str, str]:
    out: Dict[str, str] = {}
    i = start
    while i < len(tokens):
        t = tokens[i]
        if not t.startswith("-"):
            i += 1
            continue
        k = t[1:]
        if i + 1 < len(tokens) and not tokens[i + 1].startswith("-"):
            out[k] = tokens[i + 1]
            i += 2
        else:
            out[k] = "1"
            i += 1
    return out


def robust_z(x: np.ndarray) -> np.ndarray:
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    s = max(1e-9, 1.4826 * mad)
    return (x - med) / s


def safe_float(v: str, default: float = math.nan) -> float:
    try:
        return float(v)
    except Exception:
        return default


def autocorr_lag1(x: np.ndarray) -> float:
    if len(x) < 2:
        return 0.0
    a = x[:-1]
    b = x[1:]
    sa = float(np.std(a))
    sb = float(np.std(b))
    if sa <= 1e-12 or sb <= 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def switch_rate(x: np.ndarray) -> float:
    if len(x) < 3:
        return 0.0
    d = np.diff(x)
    s = np.sign(d)
    flips = np.sum((s[1:] * s[:-1]) < 0.0)
    return float(flips / max(1, len(s) - 1))


def load_gwosc_events(cache_path: Path, api_url: str) -> Dict[str, Dict]:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urlopen(api_url, timeout=30) as r:
            raw = r.read().decode("utf-8")
        cache_path.write_text(raw, encoding="utf-8")
        d = json.loads(raw)
        return d.get("events", {})
    except Exception:
        if cache_path.exists():
            d = json.loads(cache_path.read_text(encoding="utf-8"))
            return d.get("events", {})
        raise


def gps_to_mjd(gps: float) -> float:
    return float(gps) / 86400.0 + 44244.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--nanograv-archive",
        type=str,
        default="NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz",
    )
    ap.add_argument(
        "--gwosc-api",
        type=str,
        default="https://www.gw-openscience.org/eventapi/json/GWTC/",
    )
    ap.add_argument("--target-rows", type=int, default=4000)
    ap.add_argument("--heap-rows", type=int, default=12000)
    ap.add_argument("--min-snr", type=float, default=8.0)
    ap.add_argument("--max-toa-err-us", type=float, default=12.0)
    args = ap.parse_args()

    archive_path = (ROOT / args.nanograv_archive).resolve()
    if not archive_path.exists():
        raise RuntimeError(f"NANOGrav archive not found: {archive_path}")

    out_dir = ROOT / "external_confirmatory_v2" / "beta_channel_true_external"
    src_dir = out_dir / "sources"
    out_dir.mkdir(parents=True, exist_ok=True)
    src_dir.mkdir(parents=True, exist_ok=True)

    gwosc_cache = src_dir / "gwosc_gwtc_eventapi.json"
    events = load_gwosc_events(cache_path=gwosc_cache, api_url=args.gwosc_api)

    # Parse archive and keep deterministic subset by smallest hash(pair_id).
    heap: List[Tuple[int, Dict[str, object]]] = []
    n_lines_total = 0
    n_rows_seen = 0
    n_rows_quality = 0

    readme_text = ""
    with tarfile.open(archive_path, "r:gz") as tf:
        members = [m for m in tf.getmembers() if m.isfile() and "/narrowband/tim/" in m.name and m.name.endswith(".nb.tim")]
        members = sorted(members, key=lambda x: x.name)

        readme_member = next((m for m in tf.getmembers() if m.isfile() and m.name.endswith("/narrowband/README.narrowband")), None)
        if readme_member is not None:
            f = tf.extractfile(readme_member)
            if f is not None:
                readme_text = f.read().decode("utf-8", errors="ignore")
                (src_dir / "README.narrowband").write_text(readme_text, encoding="utf-8")

        for m in members:
            f = tf.extractfile(m)
            if f is None:
                continue
            pulsar = m.name.split("/")[-1].split("_PINT_")[0]
            for line_no, raw in enumerate(f, start=1):
                s = raw.decode("utf-8", errors="ignore").strip()
                if not s or s.startswith("C ") or s.startswith("FORMAT"):
                    continue
                n_lines_total += 1
                tok = s.split()
                if len(tok) < 5:
                    continue
                freq_mhz = safe_float(tok[1])
                mjd = safe_float(tok[2])
                toa_err_us = safe_float(tok[3])
                if not np.isfinite(freq_mhz) or not np.isfinite(mjd) or not np.isfinite(toa_err_us):
                    continue

                flags = parse_flag_map(tok, start=5)
                pta = flags.get("pta", "")
                if pta and pta.lower() != "nanograv":
                    continue

                snr = safe_float(flags.get("snr", "nan"))
                gof = safe_float(flags.get("gof", "nan"))
                flux = safe_float(flags.get("flux", "nan"))
                tobs = safe_float(flags.get("tobs", "nan"))
                pn = flags.get("pn", "NA")

                n_rows_seen += 1

                # Quality filter.
                if np.isfinite(snr) and snr < float(args.min_snr):
                    continue
                if toa_err_us > float(args.max_toa_err_us):
                    continue
                n_rows_quality += 1

                pair_id = f"NG15|{pulsar}|{mjd:.9f}|{freq_mhz:.3f}|{line_no}"
                h = hash_u64(pair_id)
                row = {
                    "pair_id": pair_id,
                    "pulsar": pulsar,
                    "freq_mhz": float(freq_mhz),
                    "mjd": float(mjd),
                    "toa_err_us": float(toa_err_us),
                    "snr": float(snr) if np.isfinite(snr) else math.nan,
                    "gof": float(gof) if np.isfinite(gof) else math.nan,
                    "flux": float(flux) if np.isfinite(flux) else math.nan,
                    "tobs": float(tobs) if np.isfinite(tobs) else math.nan,
                    "pn": str(pn),
                }

                # Keep smallest hashes in max-heap encoded by negative hash.
                item = (-int(h), row)
                if len(heap) < int(args.heap_rows):
                    heapq.heappush(heap, item)
                else:
                    if item[0] > heap[0][0]:
                        heapq.heapreplace(heap, item)

    if len(heap) < 1200:
        raise RuntimeError(f"Not enough quality rows after parse: {len(heap)}")

    rows = [it[1] for it in heap]
    df = pd.DataFrame(rows)

    # Deterministic sort and cut to target rows.
    df["hash"] = df["pair_id"].astype(str).map(hash_u64)
    df = df.sort_values("hash", ascending=True).head(int(args.target_rows)).copy()
    df = df.drop(columns=["hash"]).reset_index(drop=True)

    # Build theta and hxy.
    fmin = float(df["freq_mhz"].min())
    fmax = float(df["freq_mhz"].max())
    df["theta_deg"] = 180.0 * (df["freq_mhz"] - fmin) / max(1e-9, (fmax - fmin))

    x_snr = np.log1p(np.where(np.isfinite(df["snr"]), df["snr"], np.nanmedian(df["snr"])))
    x_err = np.log1p(np.clip(df["toa_err_us"].to_numpy(dtype=float), 1e-9, None))
    x_gof = np.where(np.isfinite(df["gof"]), df["gof"], np.nanmedian(df["gof"]))
    x_flux = np.log1p(np.where(np.isfinite(df["flux"]), df["flux"], np.nanmedian(df["flux"])))
    x_tobs = np.log1p(np.where(np.isfinite(df["tobs"]), df["tobs"], np.nanmedian(df["tobs"])))

    z_snr = robust_z(x_snr.astype(float))
    z_err = robust_z(x_err.astype(float))
    z_gof = robust_z(x_gof.astype(float))
    z_flux = robust_z(x_flux.astype(float))
    z_tobs = robust_z(x_tobs.astype(float))

    score = 0.45 * z_snr - 0.35 * z_err - 0.15 * np.abs(z_gof) + 0.15 * z_flux + 0.10 * z_tobs
    df["hxy"] = 1.0 / (1.0 + np.exp(-score))

    # Group scope: pulsar|pn|round(mjd,3)
    df["obs_key"] = (
        df["pulsar"].astype(str)
        + "|"
        + df["pn"].astype(str)
        + "|"
        + df["mjd"].round(3).astype(str)
    )
    df = df.sort_values(["obs_key", "theta_deg"], ascending=[True, True]).reset_index(drop=True)

    # Group-level dynamic features.
    f_std = np.zeros(len(df), dtype=float)
    f_autoc1 = np.zeros(len(df), dtype=float)
    f_switch = np.zeros(len(df), dtype=float)
    f_slope = np.zeros(len(df), dtype=float)

    for _, idx in df.groupby("obs_key").groups.items():
        ii = np.array(list(idx), dtype=int)
        y = df.loc[ii, "hxy"].to_numpy(dtype=float)
        th = df.loc[ii, "theta_deg"].to_numpy(dtype=float)
        s = float(np.std(y))
        a1 = autocorr_lag1(y)
        sw = switch_rate(y)
        if len(ii) >= 2 and float(np.std(th)) > 1e-12:
            sl = float(np.polyfit(th, y, 1)[0])
        else:
            sl = 0.0
        f_std[ii] = s
        f_autoc1[ii] = a1
        f_switch[ii] = sw
        f_slope[ii] = sl

    df["f_std"] = f_std
    df["f_autoc1"] = np.clip(f_autoc1, -1.0, 1.0)
    df["f_switch"] = np.clip(f_switch, 0.0, 1.0)
    df["f_slope"] = f_slope

    # Build intervention event table from GWOSC.
    rows_evt = []
    for eid, ed in events.items():
        gps = ed.get("GPS")
        if isinstance(gps, (int, float)):
            rows_evt.append((str(eid), float(gps)))
    if not rows_evt:
        raise RuntimeError("No GWOSC events with numeric GPS found.")
    rows_evt = sorted(rows_evt, key=lambda x: x[1])

    # Select deterministic 12 events spanning range.
    n_evt = min(12, len(rows_evt))
    idxs = np.linspace(0, len(rows_evt) - 1, num=n_evt, dtype=int)
    selected = [rows_evt[int(i)] for i in idxs]

    event_table = []
    event_mjd = {}
    for j, (eid, gps) in enumerate(selected, start=1):
        iid = f"INT_{j:02d}_{eid.replace('-', '_')}"
        mjd_evt = gps_to_mjd(gps)
        event_mjd[iid] = mjd_evt
        event_table.append(
            {
                "intervention_id": iid,
                "intervention_type": "exogenous_gw_catalog_event",
                "source_reference": "https://www.gw-openscience.org/eventapi/json/GWTC/",
                "start_utc": f"GPS:{gps:.3f}",
                "end_utc": f"GPS:{gps:.3f}",
                "is_preregistered": True,
                "is_exogenous": True,
            }
        )

    # Map each row to nearest event and regime pre/post.
    iids = list(event_mjd.keys())
    mvals = np.array([event_mjd[k] for k in iids], dtype=float)

    inter = []
    regime = []
    for mjd in df["mjd"].to_numpy(dtype=float):
        k = int(np.argmin(np.abs(mvals - mjd)))
        iid = iids[k]
        inter.append(iid)
        regime.append("pre" if mjd < mvals[k] else "post")

    df["intervention_id"] = inter
    df["regime"] = regime

    # Keep required schema first.
    out_pairs = df[
        [
            "pair_id",
            "theta_deg",
            "hxy",
            "f_std",
            "f_autoc1",
            "f_switch",
            "f_slope",
            "intervention_id",
            "regime",
        ]
    ].copy()

    n_pre = int((out_pairs["regime"].astype(str).str.lower() == "pre").sum())
    n_post = int((out_pairs["regime"].astype(str).str.lower() == "post").sum())
    if n_pre == 0 or n_post == 0:
        raise RuntimeError(f"Regime split invalid pre/post={n_pre}/{n_post}")

    beta_pairs_path = out_dir / "beta_channel_pairs.csv"
    events_path = out_dir / "intervention_events.csv"
    freeze_path = out_dir / "protocol_freeze.json"
    manifest_path = out_dir / "manifest_beta_channel.json"
    ng_copy = src_dir / archive_path.name

    out_pairs.to_csv(beta_pairs_path, index=False, quoting=csv.QUOTE_MINIMAL)
    pd.DataFrame(event_table).to_csv(events_path, index=False, quoting=csv.QUOTE_MINIMAL)

    # Keep local source copy for reproducible lineage.
    if not ng_copy.exists():
        ng_copy.write_bytes(archive_path.read_bytes())

    freeze = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "collector": "QW_1930_TRUE_EXTERNAL_BETA_CHANNEL_AUTOCOLLECTOR.py",
        "source_lineage": {
            "nanograv_archive": str(archive_path),
            "nanograv_archive_sha256": sha256_file(archive_path),
            "gwosc_api_url": args.gwosc_api,
            "gwosc_cache_sha256": sha256_file(gwosc_cache),
        },
        "protocol_lock": {
            "forbid_internal_pair_tables": True,
            "quality_filters": {
                "min_snr": float(args.min_snr),
                "max_toa_err_us": float(args.max_toa_err_us),
            },
            "row_sampling": "smallest sha256(pair_id) deterministic subset",
            "group_scope": "obs_key=pulsar|pn|round(mjd,3)",
            "intervention_mapping": "nearest GW event in MJD, regime=pre/post by signed delta",
        },
        "targets": {
            "target_rows": int(args.target_rows),
            "actual_rows": int(len(out_pairs)),
        },
    }
    freeze_path.write_text(json.dumps(freeze, ensure_ascii=False, indent=2), encoding="utf-8")

    manifest = {
        "dataset": {
            "dataset_id": f"BETA_CHANNEL_TRUE_EXTERNAL_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}",
            "provider": "NANOGRAV_GWOSC_EXTERNAL_PUBLIC",
            "externality_statement": (
                "Independent external public datasets were used: NANOGrav 15yr raw timing archive and "
                "GWOSC event catalog. Package assembled without internal analysis proxy/rebuild pair tables."
            ),
            "prepared_utc": datetime.now(timezone.utc).isoformat(),
        },
        "files": [
            {
                "role": "beta_pairs",
                "path": "beta_channel_pairs.csv",
                "sha256": sha256_file(beta_pairs_path),
            },
            {
                "role": "intervention_events",
                "path": "intervention_events.csv",
                "sha256": sha256_file(events_path),
            },
            {
                "role": "protocol_freeze",
                "path": "protocol_freeze.json",
                "sha256": sha256_file(freeze_path),
            },
            {
                "role": "source_nanograv_archive",
                "path": f"sources/{archive_path.name}",
                "sha256": sha256_file(ng_copy),
            },
            {
                "role": "source_gwosc_eventapi",
                "path": "sources/gwosc_gwtc_eventapi.json",
                "sha256": sha256_file(gwosc_cache),
            },
        ],
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "nanograv_archive": str(archive_path),
            "gwosc_api": args.gwosc_api,
        },
        "parse_stats": {
            "n_tim_lines_total": int(n_lines_total),
            "n_rows_seen": int(n_rows_seen),
            "n_rows_after_quality": int(n_rows_quality),
            "n_rows_heap_kept": int(len(heap)),
            "n_rows_output": int(len(out_pairs)),
            "n_pre": n_pre,
            "n_post": n_post,
            "n_events": int(len(event_table)),
        },
        "outputs": {
            "out_dir": str(out_dir),
            "beta_pairs": str(beta_pairs_path),
            "intervention_events": str(events_path),
            "protocol_freeze": str(freeze_path),
            "manifest": str(manifest_path),
        },
        "verdict": "TRUE_EXTERNAL_AUTOCOLLECTOR_PACKAGE_ASSEMBLED",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1930: TRUE EXTERNAL BETA-CHANNEL AUTOCOLLECTOR",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Parse Stats",
        f"- n_tim_lines_total: {out['parse_stats']['n_tim_lines_total']}",
        f"- n_rows_seen: {out['parse_stats']['n_rows_seen']}",
        f"- n_rows_after_quality: {out['parse_stats']['n_rows_after_quality']}",
        f"- n_rows_output: {out['parse_stats']['n_rows_output']}",
        f"- n_pre/n_post: {n_pre}/{n_post}",
        f"- n_events: {out['parse_stats']['n_events']}",
        "",
        "## Output Package",
        f"- {beta_pairs_path}",
        f"- {events_path}",
        f"- {freeze_path}",
        f"- {manifest_path}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1930] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1930] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1930] out_rows={len(out_pairs)} pre/post={n_pre}/{n_post}")


if __name__ == "__main__":
    main()

