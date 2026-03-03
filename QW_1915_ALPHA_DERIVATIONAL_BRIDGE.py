#!/usr/bin/env python3
"""
QW-1915: Alpha derivational bridge from constrained micromodel.

Idea:
- use constrained micromodel tuning from QW-1891 (lambda_b under nadsoliton constraints),
- map lambda_b to empirical PTA alpha through the frozen assembly scale:
    alpha = lambda_b / scale
- compare predicted alpha to empirical multisplit alpha (QW-1913).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1915_alpha_derivational_bridge.json"
OUT_MD = ROOT / "RAPORT_QW1915_ALPHA_DERIVATIONAL_BRIDGE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1891 = load("report_qw1891_derivational_constraints_from_nadsoliton.json")
    d1911 = load("report_qw1911_external_source_dataset_assembly_alpha.json")
    d1913 = load("report_qw1913_external_pta_multisplit_transfer_stress.json")

    scale = float(d1911.get("config", {}).get("scale", 0.05))
    if scale <= 0:
        raise RuntimeError("Invalid scale in QW-1911 report.")

    tuning_rows: List[Dict] = d1891.get("tuning", [])
    if not tuning_rows:
        raise RuntimeError("QW-1891 tuning rows missing.")

    lambda_rows = []
    for r in tuning_rows:
        lb = float(r["lambda_b"])
        obj = float(r["val_objective"])
        w = 1.0 / (obj + 1e-12)
        alpha_pred = lb / scale
        lambda_rows.append(
            {
                "lambda_b": lb,
                "val_objective": obj,
                "weight_inv_objective": w,
                "alpha_pred": alpha_pred,
            }
        )

    ws = np.array([r["weight_inv_objective"] for r in lambda_rows], dtype=float)
    lbs = np.array([r["lambda_b"] for r in lambda_rows], dtype=float)
    alphas = np.array([r["alpha_pred"] for r in lambda_rows], dtype=float)

    wnorm = ws / np.sum(ws)
    lambda_weighted = float(np.sum(wnorm * lbs))
    alpha_weighted = float(lambda_weighted / scale)

    selected_lambda = float(d1891.get("protocol", {}).get("selected_lambda_b", np.nan))
    alpha_selected = float(selected_lambda / scale)

    alpha_min = float(np.min(alphas))
    alpha_max = float(np.max(alphas))

    emp_summary = d1913.get("summary", {})
    alpha_emp_median = float(emp_summary.get("selected_alpha_median"))
    alpha_emp_set = [float(x) for x in (emp_summary.get("selected_alphas") or [])]
    alpha_emp_std = float(np.std(alpha_emp_set)) if alpha_emp_set else 0.0

    diff_weighted = abs(alpha_emp_median - alpha_weighted)
    diff_selected = abs(alpha_emp_median - alpha_selected)
    in_grid_range = alpha_min <= alpha_emp_median <= alpha_max

    compatibility_flags = {
        "empirical_alpha_inside_derivational_grid": bool(in_grid_range),
        "weighted_bridge_abs_diff_le_1": bool(diff_weighted <= 1.0),
        "selected_bridge_abs_diff_le_1p5": bool(diff_selected <= 1.5),
        "empirical_alpha_stable_multisplit": bool(alpha_emp_std <= 0.25),
    }

    if all(compatibility_flags.values()):
        verdict = "ALPHA_DERIVATIONAL_BRIDGE_COMPATIBLE"
    elif in_grid_range:
        verdict = "ALPHA_DERIVATIONAL_BRIDGE_PARTIAL"
    else:
        verdict = "ALPHA_DERIVATIONAL_BRIDGE_NOT_COMPATIBLE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1891_report": "report_qw1891_derivational_constraints_from_nadsoliton.json",
            "q1911_report": "report_qw1911_external_source_dataset_assembly_alpha.json",
            "q1913_report": "report_qw1913_external_pta_multisplit_transfer_stress.json",
            "scale": scale,
        },
        "derivational_mapping": {
            "formula": "alpha = lambda_b / scale",
            "tuning_rows": lambda_rows,
            "lambda_weighted_inv_objective": lambda_weighted,
            "alpha_weighted_inv_objective": alpha_weighted,
            "lambda_selected_q1891": selected_lambda,
            "alpha_selected_q1891": alpha_selected,
            "alpha_grid_range_from_q1891": [alpha_min, alpha_max],
        },
        "empirical_reference": {
            "alpha_selected_multisplit_median": alpha_emp_median,
            "alpha_selected_multisplit_values": alpha_emp_set,
            "alpha_selected_multisplit_std": alpha_emp_std,
        },
        "compatibility": {
            "abs_diff_weighted_vs_empirical": diff_weighted,
            "abs_diff_selected_vs_empirical": diff_selected,
            "flags": compatibility_flags,
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1915: ALPHA DERIVATIONAL BRIDGE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Mapping",
        f"- formula: `{out['derivational_mapping']['formula']}`",
        f"- scale: {scale}",
        f"- alpha_weighted_inv_objective: {alpha_weighted:.3f}",
        f"- alpha_selected_q1891: {alpha_selected:.3f}",
        f"- alpha_grid_range_from_q1891: [{alpha_min:.3f}, {alpha_max:.3f}]",
        "",
        "## Empirical Reference (QW-1913)",
        f"- alpha_selected_multisplit_median: {alpha_emp_median:.3f}",
        f"- alpha_selected_multisplit_values: {alpha_emp_set}",
        f"- alpha_selected_multisplit_std: {alpha_emp_std:.3f}",
        "",
        "## Compatibility",
        f"- abs_diff_weighted_vs_empirical: {diff_weighted:.3f}",
        f"- abs_diff_selected_vs_empirical: {diff_selected:.3f}",
        f"- flags: {compatibility_flags}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1915] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1915] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1915] verdict={verdict}")


if __name__ == "__main__":
    main()

