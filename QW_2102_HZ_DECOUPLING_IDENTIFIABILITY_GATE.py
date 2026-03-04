#!/usr/bin/env python3
"""
QW-2102: H(z) decoupling identifiability gate for QW-2090.

Purpose:
- check whether current H(z) input has sufficient lever arm for a strict
  two-parameter decoupling of (h0, lambda),
- keep closure claims honest by separating weak-input geometry from true
  model failure.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "h0_lambda_decoupling_input_qw2090.json"
OUT_JSON = ROOT / "report_qw2102_hz_decoupling_identifiability_gate.json"
OUT_MD = ROOT / "RAPORT_QW2102_HZ_DECOUPLING_IDENTIFIABILITY_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_metadata_complete(in_data: Dict) -> bool:
    src = str(in_data.get("source", "")).strip()
    citation = str(in_data.get("citation", "")).strip()
    ref_url = str(in_data.get("reference_url", "")).strip()
    src_sha = str(in_data.get("source_sha256", "")).strip()
    src_ver = str(in_data.get("source_version", "")).strip()
    low = f"{src} {citation} {ref_url}".lower()
    not_placeholder = bool(src) and ("placeholder" not in low)
    has_reference = bool(citation or ref_url)
    has_integrity = bool(src_sha or src_ver)
    return bool(not_placeholder and has_reference and has_integrity)


def valid_nodes(nodes: List[Dict]) -> List[Dict]:
    out: List[Dict] = []
    for row in nodes:
        try:
            z = float(row["z"])
            h = float(row["h_km_s_mpc"])
            s = abs(float(row["sigma_total"]))
        except Exception:  # noqa: BLE001
            continue
        if h <= 0.0 or s <= 0.0:
            continue
        out.append({"z": z, "h_km_s_mpc": h, "sigma_total": s})
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2102 H(z) decoupling identifiability gate")
    p.add_argument("--input", default=str(DEFAULT_INPUT), help="QW-2090 input JSON")
    args = p.parse_args()

    in_path = Path(args.input).resolve()
    input_present = in_path.exists()
    in_data = load_json(in_path) if input_present else {}

    nodes = valid_nodes(in_data.get("h_of_z_nodes", []) if input_present else [])
    omega_m = float(in_data.get("omega_m", 0.315)) if input_present else 0.315
    omega_r = float(in_data.get("omega_r", 9.0e-5)) if input_present else 9.0e-5

    z_span = None
    e_span = None
    design_condition_number = None
    e_min = None
    e_max = None
    z_min = None
    z_max = None
    fit_info = {}
    q2090_path = ROOT / "report_qw2090_h0_lambda_decoupling_gate.json"
    if q2090_path.exists():
        q2090 = load_json(q2090_path)
        fit_info = {
            "qw2090_verdict": q2090.get("verdict"),
            "qw2090_h0_rel_err_pct": q2090.get("reference", {}).get("h0_rel_err_pct"),
            "qw2090_lambda_rel_err_pct": q2090.get("reference", {}).get("lambda_rel_err_pct"),
            "qw2090_flat_h0_rel_err_pct": q2090.get("flatness_projection_diagnostic", {}).get(
                "h0_rel_err_pct"
            ),
            "qw2090_flat_lambda_rel_err_pct": q2090.get("flatness_projection_diagnostic", {}).get(
                "lambda_rel_err_pct"
            ),
        }

    if len(nodes) >= 2:
        z_vals = np.array([r["z"] for r in nodes], dtype=float)
        e_vals = np.array(
            [omega_m * (1.0 + z) ** 3 + omega_r * (1.0 + z) ** 4 for z in z_vals], dtype=float
        )
        a = np.column_stack([e_vals, np.ones_like(e_vals)])
        z_min = float(np.min(z_vals))
        z_max = float(np.max(z_vals))
        e_min = float(np.min(e_vals))
        e_max = float(np.max(e_vals))
        z_span = float(z_max - z_min)
        e_span = float(e_max - e_min)
        design_condition_number = float(np.linalg.cond(a))

    # Conservative strict-ready heuristics for two-parameter [E, 1] decoupling.
    min_nodes = 5
    min_z_span = 0.8
    min_e_span = 1.0
    max_cond = 8.0

    flags = {
        "input_present": bool(input_present),
        "source_metadata_complete": bool(source_metadata_complete(in_data) if input_present else False),
        "provenance_anchor_free": bool(in_data.get("provenance_anchor_free", False)) if input_present else False,
        "n_nodes_ge_5": bool(len(nodes) >= min_nodes),
        "z_span_ge_0p8": bool(z_span is not None and z_span >= min_z_span),
        "e_span_ge_1p0": bool(e_span is not None and e_span >= min_e_span),
        "design_condition_lt_8": bool(
            design_condition_number is not None and math.isfinite(design_condition_number) and design_condition_number < max_cond
        ),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    strict_ready = bool(all(flags.values()))
    verdict = (
        "HZ_DECOUPLING_IDENTIFIABILITY_GATE_PASS_STRICT_READY"
        if strict_ready
        else "HZ_DECOUPLING_IDENTIFIABILITY_GATE_WEAK_LEVERARM_PENDING"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "input_json": str(in_path) if input_present else None,
            "qw2090_report": "report_qw2090_h0_lambda_decoupling_gate.json" if q2090_path.exists() else None,
        },
        "thresholds": {
            "min_nodes": min_nodes,
            "min_z_span": min_z_span,
            "min_e_span": min_e_span,
            "max_design_condition_number": max_cond,
        },
        "metrics": {
            "n_nodes": len(nodes),
            "z_min": z_min,
            "z_max": z_max,
            "z_span": z_span,
            "e_min": e_min,
            "e_max": e_max,
            "e_span": e_span,
            "design_condition_number": design_condition_number,
        },
        "qw2090_context": fit_info,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "RERUN_QW2090_STRICT_DECOUPLING_WITH_EXPANDED_NONANCHOR_HZ_DATASET"
            if strict_ready
            else "EXPAND_HZ_INPUT_LEVERARM_MORE_NODES_WIDER_Z_RANGE_AND_RERUN_QW2102_QW2090"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2102: HZ DECOUPLING IDENTIFIABILITY GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Metrics",
        f"- n_nodes: {len(nodes)}",
        f"- z_span: {z_span}",
        f"- e_span: {e_span}",
        f"- design_condition_number: {design_condition_number}",
        "",
        "## Context (QW-2090)",
        f"- qw2090_verdict: {fit_info.get('qw2090_verdict')}",
        f"- qw2090_h0_rel_err_pct: {fit_info.get('qw2090_h0_rel_err_pct')}",
        f"- qw2090_lambda_rel_err_pct: {fit_info.get('qw2090_lambda_rel_err_pct')}",
        f"- qw2090_flat_h0_rel_err_pct: {fit_info.get('qw2090_flat_h0_rel_err_pct')}",
        f"- qw2090_flat_lambda_rel_err_pct: {fit_info.get('qw2090_flat_lambda_rel_err_pct')}",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2102] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2102] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2102] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

