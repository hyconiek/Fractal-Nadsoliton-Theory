#!/usr/bin/env python3
"""
QW-2106: strict external input intake gate for T3/T4.

Purpose:
- validate raw external inputs before autocollectors/gates,
- enforce explicit metadata + provenance requirements to reduce hidden
  non-strict paths.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent

DEFAULT_HZ_CSV = ROOT / "external_hz_nodes_qw2099.csv"
DEFAULT_HZ_META = ROOT / "external_hz_nodes_qw2099.metadata.json"
DEFAULT_G_JSON = ROOT / "external_gnewton_bridge_qw2101.json"
DEFAULT_G_META = ROOT / "external_gnewton_bridge_qw2101.metadata.json"

OUT_JSON = ROOT / "report_qw2106_strict_external_input_intake_gate.json"
OUT_MD = ROOT / "RAPORT_QW2106_STRICT_EXTERNAL_INPUT_INTAKE_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def metadata_complete(meta: Dict) -> bool:
    src = str(meta.get("source_dataset_name", "")).strip()
    citation = str(meta.get("citation", "")).strip()
    ref_url = str(meta.get("reference_url", "")).strip()
    ver = str(meta.get("source_version", "")).strip()
    low = f"{src} {citation} {ref_url} {ver}".lower()
    not_placeholder = bool(src and citation and ref_url and ver) and ("placeholder" not in low)
    return bool(not_placeholder)


def load_hz_rows(path: Path) -> List[Dict]:
    rows: List[Dict] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        required = {"z", "h_km_s_mpc", "sigma_total"}
        if not required.issubset(set(r.fieldnames or [])):
            raise ValueError(f"H(z) CSV must include columns: {sorted(required)}")
        for row in r:
            z = float(row["z"])
            h = float(row["h_km_s_mpc"])
            s = float(row["sigma_total"])
            rows.append({"z": z, "h": h, "s": s})
    return rows


def hz_metrics(rows: List[Dict], omega_m: float = 0.315, omega_r: float = 9e-5) -> Dict[str, float]:
    valid = [r for r in rows if r["h"] > 0.0 and r["s"] > 0.0 and math.isfinite(r["z"])]
    if not valid:
        return {
            "n_nodes_valid": 0,
            "z_span": float("nan"),
            "e_span": float("nan"),
            "design_condition_number": float("nan"),
        }
    z = np.array([r["z"] for r in valid], dtype=float)
    e = np.array([omega_m * (1.0 + zz) ** 3 + omega_r * (1.0 + zz) ** 4 for zz in z], dtype=float)
    a = np.column_stack((e, np.ones_like(e)))
    return {
        "n_nodes_valid": int(len(valid)),
        "z_span": float(np.max(z) - np.min(z)),
        "e_span": float(np.max(e) - np.min(e)),
        "design_condition_number": float(np.linalg.cond(a)),
    }


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2106 strict external input intake gate")
    p.add_argument("--hz-csv", default=str(DEFAULT_HZ_CSV), help="H(z) CSV input path.")
    p.add_argument("--hz-meta", default=str(DEFAULT_HZ_META), help="H(z) metadata JSON path.")
    p.add_argument("--g-json", default=str(DEFAULT_G_JSON), help="G bridge JSON input path.")
    p.add_argument("--g-meta", default=str(DEFAULT_G_META), help="G bridge metadata JSON path.")
    args = p.parse_args()

    hz_csv = Path(args.hz_csv).resolve()
    hz_meta_p = Path(args.hz_meta).resolve()
    g_json_p = Path(args.g_json).resolve()
    g_meta_p = Path(args.g_meta).resolve()

    hz_csv_exists = hz_csv.exists()
    hz_meta_exists = hz_meta_p.exists()
    g_json_exists = g_json_p.exists()
    g_meta_exists = g_meta_p.exists()

    hz_rows: List[Dict] = []
    hz_meta: Dict = {}
    g_data: Dict = {}
    g_meta: Dict = {}
    hz_error = ""
    g_error = ""

    if hz_csv_exists:
        try:
            hz_rows = load_hz_rows(hz_csv)
        except Exception as exc:  # noqa: BLE001
            hz_error = str(exc)
    if hz_meta_exists:
        try:
            hz_meta = load_json(hz_meta_p)
        except Exception as exc:  # noqa: BLE001
            hz_error = hz_error or str(exc)
    if g_json_exists:
        try:
            g_data = load_json(g_json_p)
        except Exception as exc:  # noqa: BLE001
            g_error = str(exc)
    if g_meta_exists:
        try:
            g_meta = load_json(g_meta_p)
        except Exception as exc:  # noqa: BLE001
            g_error = g_error or str(exc)

    m = hz_metrics(hz_rows)

    hz_flags = {
        "hz_csv_exists": hz_csv_exists,
        "hz_csv_parse_ok": bool(hz_csv_exists and not hz_error),
        "hz_meta_exists": hz_meta_exists,
        "hz_meta_complete": bool(hz_meta_exists and metadata_complete(hz_meta)),
        "hz_provenance_anchor_free": bool(hz_meta.get("provenance_anchor_free", False)),
        "hz_n_nodes_ge_5": bool(m["n_nodes_valid"] >= 5),
        "hz_z_span_ge_0p8": bool(math.isfinite(m["z_span"]) and m["z_span"] >= 0.8),
        "hz_e_span_ge_1p0": bool(math.isfinite(m["e_span"]) and m["e_span"] >= 1.0),
        "hz_design_condition_lt_8": bool(
            math.isfinite(m["design_condition_number"]) and m["design_condition_number"] < 8.0
        ),
    }

    g_origin = str(g_data.get("bridge_observable_origin", g_meta.get("bridge_observable_origin", ""))).strip()
    g_dimless = g_data.get("g_dimensionless_mu_ref", None)
    g_si = g_data.get("g_si", None)
    g_flags = {
        "g_json_exists": g_json_exists,
        "g_json_parse_ok": bool(g_json_exists and not g_error),
        "g_meta_exists": g_meta_exists,
        "g_meta_complete": bool(g_meta_exists and metadata_complete(g_meta)),
        "g_provenance_anchor_free": bool(g_meta.get("provenance_anchor_free", False)),
        "g_not_seeded_from_registry": bool(not g_meta.get("seeded_from_registry", False)),
        "g_origin_external_dimensionless": bool(g_origin == "external_dimensionless_observable"),
        "g_dimensionless_present_positive": bool(
            isinstance(g_dimless, (int, float)) and math.isfinite(float(g_dimless)) and float(g_dimless) > 0.0
        ),
        "g_si_not_primary": bool(g_si is None),
    }

    hz_ready = all(hz_flags.values())
    g_ready = all(g_flags.values())
    all_ready = bool(hz_ready and g_ready)

    if all_ready:
        verdict = "STRICT_EXTERNAL_INPUT_INTAKE_GATE_PASS"
    else:
        verdict = "STRICT_EXTERNAL_INPUT_INTAKE_GATE_PENDING"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "hz_csv": str(hz_csv),
            "hz_meta": str(hz_meta_p),
            "g_json": str(g_json_p),
            "g_meta": str(g_meta_p),
        },
        "errors": {
            "hz_error": hz_error,
            "g_error": g_error,
        },
        "hz_metrics": m,
        "hz_flags": hz_flags,
        "g_flags": g_flags,
        "pass_count": int(sum(1 for v in {**hz_flags, **g_flags}.values() if bool(v))),
        "total_flags": int(len(hz_flags) + len(g_flags)),
        "verdict": verdict,
        "required_next_step": (
            "RUN_QW2099_QW2101_QW2102_QW2103_QW2090_QW2092_QW2104_QW2105"
            if all_ready
            else "FILL_METADATA_AND_STRICT_EXTERNAL_INPUT_GAPS_THEN_RERUN_QW2106"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2106: STRICT EXTERNAL INPUT INTAKE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {out['pass_count']}/{out['total_flags']}",
        "",
        "## H(z) Metrics",
        f"- n_nodes_valid: {m['n_nodes_valid']}",
        f"- z_span: {m['z_span']}",
        f"- e_span: {m['e_span']}",
        f"- design_condition_number: {m['design_condition_number']}",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2106] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2106] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2106] verdict={verdict} pass_count={out['pass_count']}/{out['total_flags']}")


if __name__ == "__main__":
    main()
