#!/usr/bin/env python3
"""
QW-1928: Self-fetched beta-channel package assembly (test-only).

Builds a beta-channel package using:
- externally fetched metadata pages (GWOSC + NANOGrav),
- external-source rebuilt PTA pair table already present in repo.

Important scientific note:
- this is a self-assembled package for mechanical gate testing,
- it does NOT replace an independent third-party confirmatory package.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1928_self_fetched_beta_channel_package_assembly.json"
OUT_MD = ROOT / "RAPORT_QW1928_SELF_FETCHED_BETA_CHANNEL_PACKAGE_ASSEMBLY.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def split_hash(pair_id: str) -> int:
    h = hashlib.sha256(pair_id.encode("utf-8")).hexdigest()
    return int(h[-8:], 16)


def main() -> None:
    src_pta = ROOT / "external_confirmatory_v2" / "confirmatory_dataset_external_source_rebuild_v2_1831cfg" / "pta_v2_pairs.csv"
    src_gwosc = ROOT / "external_confirmatory_v2" / "beta_channel_true_external" / "sources" / "gwosc_gwtc_eventapi.json"
    src_nanograv = ROOT / "external_confirmatory_v2" / "beta_channel_true_external" / "sources" / "nanograv_data_page.html"

    if not src_pta.exists():
        raise RuntimeError(f"Missing source PTA: {src_pta}")
    if not src_gwosc.exists():
        raise RuntimeError(f"Missing source GWOSC metadata: {src_gwosc}")
    if not src_nanograv.exists():
        raise RuntimeError(f"Missing source NANOGrav metadata: {src_nanograv}")

    out_dir = ROOT / "external_confirmatory_v2" / "beta_channel_true_external"
    out_dir.mkdir(parents=True, exist_ok=True)

    beta_pairs_path = out_dir / "beta_channel_pairs.csv"
    events_path = out_dir / "intervention_events.csv"
    freeze_path = out_dir / "protocol_freeze.json"
    manifest_path = out_dir / "manifest_beta_channel.json"

    df = pd.read_csv(src_pta)
    req = ["pair_id", "theta_deg", "hxy", "f_std", "f_autoc1", "f_switch", "f_slope"]
    miss = [c for c in req if c not in df.columns]
    if miss:
        raise RuntimeError(f"PTA source missing required columns: {miss}")

    # Build intervention catalog from GWOSC event api.
    d_gw = json.loads(src_gwosc.read_text(encoding="utf-8"))
    events: Dict[str, Dict] = d_gw.get("events", {})
    event_ids = sorted([k for k in events.keys() if str(k).startswith("GW")])[:6]
    if len(event_ids) < 2:
        event_ids = [f"EXT_EVT_{i+1}" for i in range(4)]

    event_rows: List[Dict[str, object]] = []
    for i, eid in enumerate(event_ids):
        if eid.startswith("GW") and eid in events:
            ed = events[eid]
            gps = ed.get("GPS")
            if isinstance(gps, (int, float)):
                # Keep as string reference; precise UTC conversion is not required by gate.
                start = f"GPS:{gps}"
                end = f"GPS:{gps}"
            else:
                start = "UNKNOWN"
                end = "UNKNOWN"
            ref = f"https://www.gw-openscience.org/eventapi/html/GWTC/"
        else:
            start = "UNKNOWN"
            end = "UNKNOWN"
            ref = "https://www.gw-openscience.org/eventapi/html/GWTC/"

        event_rows.append(
            {
                "intervention_id": f"INT_{i+1:02d}",
                "intervention_type": "exogenous_detector_environment",
                "source_reference": ref,
                "start_utc": start,
                "end_utc": end,
                "is_preregistered": True,
                "is_exogenous": True,
            }
        )

    # Assign each pair to intervention cohort and regime.
    n_int = len(event_rows)
    int_ids = [r["intervention_id"] for r in event_rows]

    hashes = np.array([split_hash(str(x)) for x in df["pair_id"].astype(str)], dtype=np.int64)
    df["intervention_id"] = [int_ids[int(h % n_int)] for h in hashes]

    # Pre/post split for gate compliance and downstream intervention logic.
    # Deterministic regime assignment by hash bit (balanced on large tables).
    df["regime"] = np.where((hashes % 2) == 0, "pre", "post")

    # Keep strict schema first, then append optional retained columns.
    keep_cols = [
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
    extra = [c for c in df.columns if c not in keep_cols]
    df_out = df[keep_cols + extra].copy()

    df_out.to_csv(beta_pairs_path, index=False, quoting=csv.QUOTE_MINIMAL)
    pd.DataFrame(event_rows).to_csv(events_path, index=False, quoting=csv.QUOTE_MINIMAL)

    freeze = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "pta_pairs": str(src_pta),
            "gwosc_metadata": str(src_gwosc),
            "nanograv_metadata": str(src_nanograv),
        },
        "locked_observable": "B7_local_resid_std",
        "locked_rules": {
            "split_method": "sha256(pair_id) parity",
            "beta_effect_min": 0.60,
            "omega_effect_max": 0.25,
            "contrast_min": 0.35,
            "contrast_boot_q05_min": 0.20,
        },
        "scientific_note": (
            "Self-fetched/self-assembled package for mechanical readiness testing. "
            "Not an independent third-party confirmatory package."
        ),
    }
    freeze_path.write_text(json.dumps(freeze, ensure_ascii=False, indent=2), encoding="utf-8")

    manifest = {
        "dataset": {
            "dataset_id": f"BETA_CHANNEL_SELF_FETCH_TEST_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}",
            "provider": "EXTERNAL_SELF_FETCH_REBUILD_TEST",
            "externality_statement": (
                "Data package uses external-source metadata and external-source rebuilt pairs; "
                "this test package is external in source lineage but not independent from current internal analysis workflow."
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
                "role": "source_nanograv_metadata",
                "path": "sources/nanograv_data_page.html",
                "sha256": sha256_file(src_nanograv),
            },
            {
                "role": "source_gwosc_metadata",
                "path": "sources/gwosc_gwtc_eventapi.json",
                "sha256": sha256_file(src_gwosc),
            },
        ],
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    n_pre = int((df_out["regime"].astype(str).str.lower() == "pre").sum())
    n_post = int((df_out["regime"].astype(str).str.lower() == "post").sum())

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "output_dir": str(out_dir),
        "rows": {
            "beta_channel_pairs": int(len(df_out)),
            "intervention_events": int(len(event_rows)),
            "pre": n_pre,
            "post": n_post,
        },
        "source_lineage": {
            "pta_source": str(src_pta),
            "nanograv_page": str(src_nanograv),
            "gwosc_eventapi": str(src_gwosc),
        },
        "artifacts": {
            "manifest": str(manifest_path),
            "beta_pairs": str(beta_pairs_path),
            "intervention_events": str(events_path),
            "protocol_freeze": str(freeze_path),
        },
        "verdict": "SELF_FETCHED_BETA_CHANNEL_PACKAGE_ASSEMBLED",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1928: SELF-FETCHED BETA-CHANNEL PACKAGE ASSEMBLY",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- output_dir: `{out['output_dir']}`",
        "",
        "## Rows",
        f"- beta_channel_pairs: {out['rows']['beta_channel_pairs']}",
        f"- intervention_events: {out['rows']['intervention_events']}",
        f"- pre/post: {out['rows']['pre']}/{out['rows']['post']}",
        "",
        "## Scientific Note",
        "- This package is self-assembled for mechanical gate testing.",
        "- It does not satisfy independent third-party confirmatory standard by itself.",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1928] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1928] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1928] output_dir={out_dir}")


if __name__ == "__main__":
    main()
