#!/usr/bin/env python3
"""
QW-2107: deterministic strict-design search for external H(z) node expansion.

Purpose:
- do not fabricate measurements,
- compute a minimal redshift acquisition design that can satisfy strict
  identifiability thresholds used by QW-2099/QW-2102/QW-2106,
- produce an auditable acquisition target list for external teams.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_HZ_CSV = ROOT / "external_hz_nodes_qw2099.csv"
OUT_JSON = ROOT / "report_qw2107_hz_strict_design_search.json"
OUT_MD = ROOT / "RAPORT_QW2107_HZ_STRICT_DESIGN_SEARCH.md"


def load_z_nodes(path: Path) -> List[float]:
    vals: List[float] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        required = {"z", "h_km_s_mpc", "sigma_total"}
        if not required.issubset(set(r.fieldnames or [])):
            raise ValueError(f"CSV must contain columns: {sorted(required)}")
        for row in r:
            z = float(row["z"])
            h = float(row["h_km_s_mpc"])
            s = float(row["sigma_total"])
            if h > 0.0 and s > 0.0 and math.isfinite(z):
                vals.append(float(z))
    vals = sorted(set(vals))
    if len(vals) < 3:
        raise ValueError("Need at least 3 valid H(z) nodes for design search.")
    return vals


def metrics_for_z(
    z_values: Sequence[float],
    omega_m: float,
    omega_r: float,
) -> Dict[str, float]:
    z = np.array(sorted(z_values), dtype=float)
    e = np.array([omega_m * (1.0 + zz) ** 3 + omega_r * (1.0 + zz) ** 4 for zz in z], dtype=float)
    a = np.column_stack((e, np.ones_like(e)))
    return {
        "n_nodes": int(len(z)),
        "z_min": float(np.min(z)),
        "z_max": float(np.max(z)),
        "z_span": float(np.max(z) - np.min(z)),
        "e_min": float(np.min(e)),
        "e_max": float(np.max(e)),
        "e_span": float(np.max(e) - np.min(e)),
        "design_condition_number": float(np.linalg.cond(a)),
    }


def strict_flags(m: Dict[str, float], thr: Dict[str, float]) -> Dict[str, bool]:
    return {
        "n_nodes_ge_min": bool(m["n_nodes"] >= int(thr["min_nodes"])),
        "z_span_ge_min": bool(m["z_span"] >= float(thr["min_z_span"])),
        "e_span_ge_min": bool(m["e_span"] >= float(thr["min_e_span"])),
        "design_condition_lt_max": bool(
            math.isfinite(m["design_condition_number"])
            and m["design_condition_number"] < float(thr["max_design_condition_number"])
        ),
    }


def all_pass(flags: Dict[str, bool]) -> bool:
    return all(bool(v) for v in flags.values())


def build_candidate_grid(
    z_min: float,
    z_max: float,
    step: float,
    existing: Sequence[float],
    min_sep: float,
) -> List[float]:
    out: List[float] = []
    n = int(round((z_max - z_min) / step))
    for i in range(n + 1):
        z = round(z_min + i * step, 8)
        if any(abs(z - e) < min_sep for e in existing):
            continue
        out.append(z)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2107 H(z) strict design search")
    p.add_argument("--hz-csv", default=str(DEFAULT_HZ_CSV), help="Input H(z) CSV")
    p.add_argument("--omega-m", type=float, default=0.315)
    p.add_argument("--omega-r", type=float, default=9.0e-5)

    p.add_argument("--min-nodes", type=int, default=5)
    p.add_argument("--min-z-span", type=float, default=0.8)
    p.add_argument("--min-e-span", type=float, default=1.0)
    p.add_argument("--max-cond", type=float, default=8.0)

    p.add_argument("--candidate-z-min", type=float, default=0.1)
    p.add_argument("--candidate-z-max", type=float, default=3.0)
    p.add_argument("--candidate-step", type=float, default=0.02)
    p.add_argument("--min-separation", type=float, default=0.015)
    p.add_argument("--max-added-nodes", type=int, default=4)
    p.add_argument("--top-k", type=int, default=20)
    args = p.parse_args()

    hz_csv = Path(args.hz_csv).resolve()
    if not hz_csv.exists():
        raise FileNotFoundError(f"H(z) CSV not found: {hz_csv}")

    existing = load_z_nodes(hz_csv)
    thresholds = {
        "min_nodes": int(args.min_nodes),
        "min_z_span": float(args.min_z_span),
        "min_e_span": float(args.min_e_span),
        "max_design_condition_number": float(args.max_cond),
    }

    base_metrics = metrics_for_z(existing, omega_m=float(args.omega_m), omega_r=float(args.omega_r))
    base_flags = strict_flags(base_metrics, thresholds)

    k_min = max(0, int(args.min_nodes) - len(existing))
    k_max = max(k_min, int(args.max_added_nodes))

    candidates = build_candidate_grid(
        z_min=float(args.candidate_z_min),
        z_max=float(args.candidate_z_max),
        step=float(args.candidate_step),
        existing=existing,
        min_sep=float(args.min_separation),
    )

    found: List[Dict] = []
    search_summary: List[Dict] = []

    for k in range(k_min, k_max + 1):
        n_tested = 0
        n_pass = 0
        best_cond = None
        for combo in itertools.combinations(candidates, k):
            n_tested += 1
            z_all = sorted(existing + list(combo))
            m = metrics_for_z(z_all, omega_m=float(args.omega_m), omega_r=float(args.omega_r))
            f = strict_flags(m, thresholds)
            if best_cond is None or m["design_condition_number"] < best_cond:
                best_cond = m["design_condition_number"]
            if all_pass(f):
                n_pass += 1
                found.append(
                    {
                        "added_nodes": [float(x) for x in combo],
                        "metrics": m,
                        "flags": f,
                    }
                )
        search_summary.append(
            {
                "k_added": int(k),
                "n_tested": int(n_tested),
                "n_pass": int(n_pass),
                "best_condition_number_seen": None if best_cond is None else float(best_cond),
            }
        )
        if n_pass > 0:
            break

    def rank_key(item: Dict) -> tuple:
        added = item["added_nodes"]
        m = item["metrics"]
        return (
            len(added),
            max(added) if added else -1.0,
            m["design_condition_number"],
            m["z_span"],
            m["e_span"],
        )

    found_sorted = sorted(found, key=rank_key)
    top = found_sorted[: int(args.top_k)]

    if top:
        verdict = "HZ_STRICT_DESIGN_FOUND"
        required_next_step = "COLLECT_REAL_HZ_MEASUREMENTS_AT_RECOMMENDED_Z_AND_RERUN_QW2106_QW2099_QW2102_QW2090"
    else:
        verdict = "HZ_STRICT_DESIGN_NOT_FOUND_IN_SEARCH_GRID"
        required_next_step = "EXPAND_SEARCH_GRID_OR_RELAX_NONPHYSICAL_CONSTRAINTS_ONLY_WITH_EXPLICIT_GOVERNANCE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_hz_csv": str(hz_csv),
        "existing_z_nodes": existing,
        "omega_m": float(args.omega_m),
        "omega_r": float(args.omega_r),
        "thresholds": thresholds,
        "base_metrics": base_metrics,
        "base_flags": base_flags,
        "candidate_grid": {
            "z_min": float(args.candidate_z_min),
            "z_max": float(args.candidate_z_max),
            "step": float(args.candidate_step),
            "min_separation": float(args.min_separation),
            "n_candidates": int(len(candidates)),
        },
        "k_min_required_by_node_count": int(k_min),
        "k_max_searched": int(k_max),
        "search_summary": search_summary,
        "n_solutions": int(len(found_sorted)),
        "top_solutions": top,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2107: HZ STRICT DESIGN SEARCH",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- input_hz_csv: `{hz_csv}`",
        f"- existing_nodes: `{len(existing)}`",
        f"- threshold_min_nodes: `{thresholds['min_nodes']}`",
        f"- threshold_min_z_span: `{thresholds['min_z_span']}`",
        f"- threshold_min_e_span: `{thresholds['min_e_span']}`",
        f"- threshold_max_cond: `{thresholds['max_design_condition_number']}`",
        "",
        "## Base Metrics",
        f"- n_nodes: `{base_metrics['n_nodes']}`",
        f"- z_span: `{base_metrics['z_span']:.6g}`",
        f"- e_span: `{base_metrics['e_span']:.6g}`",
        f"- design_condition_number: `{base_metrics['design_condition_number']:.6g}`",
        "",
        "## Search",
        f"- candidate_grid_size: `{len(candidates)}`",
        f"- k_min_required: `{k_min}`",
        f"- k_max_searched: `{k_max}`",
        f"- n_solutions: `{len(found_sorted)}`",
        "",
        "## Top Solutions (z additions)",
    ]
    if top:
        for i, item in enumerate(top[:10], start=1):
            add = ", ".join(f"{z:.3f}" for z in item["added_nodes"])
            m = item["metrics"]
            lines.append(
                f"{i}. add [{add}] -> z_span={m['z_span']:.6g}, e_span={m['e_span']:.6g}, cond={m['design_condition_number']:.6g}"
            )
    else:
        lines.append("- no strict solution found in this search grid")

    lines += ["", "## Artifact", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2107] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2107] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2107] verdict={verdict} solutions={len(found_sorted)}")


if __name__ == "__main__":
    main()
