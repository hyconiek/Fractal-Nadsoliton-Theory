#!/usr/bin/env python3
"""
QW-2099: external H(z) decoupling autocollector for QW-2090.

Builds a strict-input artifact with provenance metadata:
- h0_lambda_decoupling_input_qw2090.json
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_CSV = ROOT / "external_hz_nodes_qw2099.csv"
DEFAULT_OUT_INPUT = ROOT / "h0_lambda_decoupling_input_qw2090.json"
OUT_JSON = ROOT / "report_qw2099_hz_external_decoupling_autocollector.json"
OUT_MD = ROOT / "RAPORT_QW2099_HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTOR.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_nodes(csv_path: Path) -> List[Dict]:
    rows: List[Dict] = []
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        need = {"z", "h_km_s_mpc", "sigma_total"}
        if not need.issubset(set(r.fieldnames or [])):
            raise ValueError(f"CSV must contain columns: {sorted(need)}")
        for row in r:
            z = float(row["z"])
            h = float(row["h_km_s_mpc"])
            s = float(row["sigma_total"])
            if h <= 0.0 or s <= 0.0:
                continue
            rows.append({"z": z, "h_km_s_mpc": h, "sigma_total": s})
    if len(rows) < 3:
        raise ValueError("Need at least 3 valid H(z) nodes.")
    return rows


def identifiability_metrics(nodes: List[Dict], omega_m: float, omega_r: float) -> Dict[str, float]:
    z_vals = np.array([float(r["z"]) for r in nodes], dtype=float)
    e_vals = np.array(
        [omega_m * (1.0 + float(r["z"])) ** 3 + omega_r * (1.0 + float(r["z"])) ** 4 for r in nodes],
        dtype=float,
    )
    a = np.column_stack((e_vals, np.ones_like(e_vals)))
    cond = float(np.linalg.cond(a))
    return {
        "n_nodes": int(len(nodes)),
        "z_min": float(np.min(z_vals)),
        "z_max": float(np.max(z_vals)),
        "z_span": float(np.max(z_vals) - np.min(z_vals)),
        "e_min": float(np.min(e_vals)),
        "e_max": float(np.max(e_vals)),
        "e_span": float(np.max(e_vals) - np.min(e_vals)),
        "design_condition_number": cond,
    }


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2099 external H(z) decoupling autocollector")
    p.add_argument("--nodes-csv", default=str(DEFAULT_CSV), help="CSV with columns z,h_km_s_mpc,sigma_total")
    p.add_argument("--output-input", default=str(DEFAULT_OUT_INPUT), help="Output JSON for QW-2090")
    p.add_argument("--source", default="external_hz_nodes_autocollected", help="Short source label")
    p.add_argument("--citation", required=True, help="Human-readable source citation")
    p.add_argument("--reference-url", required=True, help="Reference URL or DOI resolver")
    p.add_argument("--source-version", required=True, help="Dataset/paper version or release tag")
    p.add_argument("--omega-m", type=float, default=0.315, help="Omega_m used in decoupling model")
    p.add_argument("--omega-r", type=float, default=9.0e-5, help="Omega_r used in decoupling model")
    p.add_argument(
        "--require-strict-ready",
        action="store_true",
        help="Return non-zero exit code unless strict H(z) identifiability thresholds are met.",
    )
    args = p.parse_args()

    csv_path = Path(args.nodes_csv).resolve()
    out_input = Path(args.output_input).resolve()
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    nodes = load_nodes(csv_path)
    sha = sha256_file(csv_path)
    omega_m = float(args.omega_m)
    omega_r = float(args.omega_r)

    metrics = identifiability_metrics(nodes, omega_m=omega_m, omega_r=omega_r)
    thresholds = {
        "min_nodes": 5,
        "min_z_span": 0.8,
        "min_e_span": 1.0,
        "max_design_condition_number": 8.0,
    }
    strict_flags = {
        "n_nodes_ge_5": bool(metrics["n_nodes"] >= thresholds["min_nodes"]),
        "z_span_ge_0p8": bool(metrics["z_span"] >= thresholds["min_z_span"]),
        "e_span_ge_1p0": bool(metrics["e_span"] >= thresholds["min_e_span"]),
        "design_condition_lt_8": bool(
            math.isfinite(metrics["design_condition_number"])
            and metrics["design_condition_number"] < thresholds["max_design_condition_number"]
        ),
    }
    strict_ready = bool(all(strict_flags.values()))

    payload = {
        "source": str(args.source),
        "citation": str(args.citation),
        "reference_url": str(args.reference_url),
        "source_version": str(args.source_version),
        "source_sha256": sha,
        "provenance_anchor_free": True,
        "omega_m": omega_m,
        "omega_r": omega_r,
        "h_of_z_nodes": nodes,
    }
    out_input.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    z_vals = [r["z"] for r in nodes]
    h_vals = [r["h_km_s_mpc"] for r in nodes]
    verdict = (
        "HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTED_STRICT_READY"
        if strict_ready
        else "HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTED_WEAK_LEVERARM"
    )
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": verdict,
        "nodes_csv": str(csv_path),
        "nodes_csv_sha256": sha,
        "n_nodes": len(nodes),
        "z_range": [float(min(z_vals)), float(max(z_vals))],
        "h_range_km_s_mpc": [float(min(h_vals)), float(max(h_vals))],
        "strict_thresholds": thresholds,
        "strict_identifiability_metrics": metrics,
        "strict_identifiability_flags": strict_flags,
        "strict_ready": strict_ready,
        "output_input_json": str(out_input),
        "required_next_step": (
            f"RUN_QW2102_AND_QW2090_WITH:{out_input.name}"
            if strict_ready
            else "EXPAND_HZ_NODES_TO_STRICT_THRESHOLD_AND_RERUN_QW2099"
        ),
    }
    OUT_JSON.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2099: HZ EXTERNAL DECOUPLING AUTOCOLLECTOR",
        "",
        f"- Date UTC: {report['generated_utc']}",
        f"- Verdict: **{report['verdict']}**",
        f"- nodes_csv: `{csv_path}`",
        f"- nodes_csv_sha256: `{sha}`",
        f"- n_nodes: `{len(nodes)}`",
        f"- z_range: `{report['z_range'][0]:.6g} .. {report['z_range'][1]:.6g}`",
        f"- h_range_km_s_mpc: `{report['h_range_km_s_mpc'][0]:.6g} .. {report['h_range_km_s_mpc'][1]:.6g}`",
        f"- strict_ready: `{strict_ready}`",
        "",
        "## Strict Identifiability",
        f"- n_nodes_ge_5: `{strict_flags['n_nodes_ge_5']}` (n={metrics['n_nodes']})",
        f"- z_span_ge_0p8: `{strict_flags['z_span_ge_0p8']}` (z_span={metrics['z_span']:.6g})",
        f"- e_span_ge_1p0: `{strict_flags['e_span_ge_1p0']}` (e_span={metrics['e_span']:.6g})",
        f"- design_condition_lt_8: `{strict_flags['design_condition_lt_8']}` "
        f"(cond={metrics['design_condition_number']:.6g})",
        "",
        "## Output",
        f"- input JSON: `{out_input}`",
        f"- report JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2099] Saved input JSON: {out_input}")
    print(f"[QW-2099] Saved report JSON: {OUT_JSON.name}")
    print(f"[QW-2099] Saved report MD:   {OUT_MD.name}")
    print(f"[QW-2099] verdict={report['verdict']} n_nodes={len(nodes)} strict_ready={strict_ready}")

    if args.require_strict_ready and not strict_ready:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
