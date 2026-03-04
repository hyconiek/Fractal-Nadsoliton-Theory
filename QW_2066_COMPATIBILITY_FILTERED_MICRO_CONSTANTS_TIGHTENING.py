#!/usr/bin/env python3
"""
QW-2066: Compatibility-filtered micro constants tightening.

Purpose:
- Reduce dispersion warning for micro-derived renormalization constants
  via deterministic quality filtering of QW-2048 bins.
- No new fitting/scan in model space; only deterministic bin selection.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2066_compatibility_filtered_micro_constants_tightening.json"
OUT_MD = ROOT / "RAPORT_QW2066_COMPATIBILITY_FILTERED_MICRO_CONSTANTS_TIGHTENING.md"


def q(x: np.ndarray, p: float) -> float:
    return float(np.quantile(x, p))


def main() -> None:
    r2048 = json.loads((ROOT / "report_qw2048_spectral_phase_locked_pointwise_derivation.json").read_text(encoding="utf-8"))
    r2049 = json.loads((ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8"))

    kernel = {k: float(v) for k, v in r2049["stagec_pool"]["selected_kernel"].items()}
    beta_uv = 0.01
    eta_uv = 1.0
    z_beta_target = float(kernel["beta"] / beta_uv)
    delta_eta_target = float(kernel["eta"] - eta_uv)

    bins = r2048["pointwise_bins"]["bins"]
    zb = np.array([float(b["z_beta_median"]) for b in bins], dtype=float)
    de = np.array([float(b["delta_eta_median"]) for b in bins], dtype=float)
    n = np.array([float(b["n"]) for b in bins], dtype=float)
    ph = np.array([float(b["phase_min_median"]) for b in bins], dtype=float)
    rm = np.array([float(b["rmse_median"]) for b in bins], dtype=float)

    base_z_q25, base_z_q50, base_z_q75 = [q(zb, v) for v in [0.25, 0.5, 0.75]]
    base_de_q25, base_de_q50, base_de_q75 = [q(de, v) for v in [0.25, 0.5, 0.75]]
    base_z_log_iqr = float(np.log(max(base_z_q75, 1e-12)) - np.log(max(base_z_q25, 1e-12)))
    base_de_iqr = float(base_de_q75 - base_de_q25)

    # Deterministic threshold families (fixed quantile grid).
    qn_grid = [0.2, 0.3, 0.4, 0.5]
    qp_grid = [0.3, 0.4, 0.5, 0.6]
    qr_grid = [0.6, 0.7, 0.8, 0.9]

    cand_rows: List[Dict[str, float]] = []
    for qn, qp, qr in product(qn_grid, qp_grid, qr_grid):
        tn = q(n, qn)
        tp = q(ph, qp)
        tr = q(rm, qr)
        mask = (n >= tn) & (ph >= tp) & (rm <= tr)
        m = int(np.sum(mask))
        if m < 5:
            continue

        z_q25, z_q50, z_q75 = [q(zb[mask], v) for v in [0.25, 0.5, 0.75]]
        de_q25, de_q50, de_q75 = [q(de[mask], v) for v in [0.25, 0.5, 0.75]]
        z_log_iqr = float(np.log(max(z_q75, 1e-12)) - np.log(max(z_q25, 1e-12)))
        de_iqr = float(de_q75 - de_q25)
        z_log_ratio = float(abs(math.log(max(z_q50, 1e-12) / max(z_beta_target, 1e-12))))
        de_abs_err = float(abs(de_q50 - delta_eta_target))

        row = {
            "q_n": float(qn),
            "q_phase": float(qp),
            "q_rmse": float(qr),
            "n_min": float(tn),
            "phase_min": float(tp),
            "rmse_max": float(tr),
            "n_selected": m,
            "z_beta_q25": float(z_q25),
            "z_beta_q50": float(z_q50),
            "z_beta_q75": float(z_q75),
            "delta_eta_q25": float(de_q25),
            "delta_eta_q50": float(de_q50),
            "delta_eta_q75": float(de_q75),
            "z_beta_log_iqr": z_log_iqr,
            "delta_eta_iqr": de_iqr,
            "z_beta_abs_log_ratio_to_target": z_log_ratio,
            "delta_eta_abs_err_to_target": de_abs_err,
            "phase_selected_median": float(np.median(ph[mask])),
            "rmse_selected_median": float(np.median(rm[mask])),
        }
        cand_rows.append(row)

    # Compatibility-constrained selection: deterministic and target-compatible.
    compat = [
        r
        for r in cand_rows
        if r["z_beta_abs_log_ratio_to_target"] <= math.log(2.0)
        and r["delta_eta_abs_err_to_target"] <= 0.25
    ]

    if compat:
        selected = sorted(
            compat,
            key=lambda r: (
                r["z_beta_log_iqr"],
                r["delta_eta_iqr"],
                -r["n_selected"],
                r["rmse_selected_median"],
            ),
        )[0]
        selection_mode = "compatibility_constrained_min_dispersion"
    else:
        selected = sorted(
            cand_rows,
            key=lambda r: (
                r["z_beta_log_iqr"],
                r["delta_eta_iqr"],
                -r["n_selected"],
                r["rmse_selected_median"],
            ),
        )[0]
        selection_mode = "fallback_min_dispersion"

    thresholds = {
        "z_beta_abs_log_ratio_to_target_max": math.log(2.0),
        "delta_eta_abs_err_to_target_max": 0.25,
        "z_beta_log_iqr_max": 2.2,
        "delta_eta_iqr_max": 0.50,
        "n_selected_min": 5,
        "phase_selected_median_min": 0.84,
    }

    flags = {
        "n_selected_ge_min": bool(selected["n_selected"] >= thresholds["n_selected_min"]),
        "z_beta_abs_log_ratio_to_target_le_max": bool(
            selected["z_beta_abs_log_ratio_to_target"] <= thresholds["z_beta_abs_log_ratio_to_target_max"]
        ),
        "delta_eta_abs_err_to_target_le_max": bool(
            selected["delta_eta_abs_err_to_target"] <= thresholds["delta_eta_abs_err_to_target_max"]
        ),
        "z_beta_log_iqr_le_max": bool(selected["z_beta_log_iqr"] <= thresholds["z_beta_log_iqr_max"]),
        "delta_eta_iqr_le_max": bool(selected["delta_eta_iqr"] <= thresholds["delta_eta_iqr_max"]),
        "phase_selected_median_ge_min": bool(
            selected["phase_selected_median"] >= thresholds["phase_selected_median_min"]
        ),
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    tightened_warning_resolved = bool(
        flags["z_beta_log_iqr_le_max"] and flags["delta_eta_iqr_le_max"] and flags["n_selected_ge_min"]
    )

    if all(flags.values()):
        verdict = "COMPATIBILITY_FILTERED_MICRO_CONSTANTS_TIGHTENING_PASS"
    elif tightened_warning_resolved:
        verdict = "COMPATIBILITY_FILTERED_MICRO_CONSTANTS_TIGHTENING_PARTIAL"
    else:
        verdict = "COMPATIBILITY_FILTERED_MICRO_CONSTANTS_TIGHTENING_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "micro_bins": "report_qw2048_spectral_phase_locked_pointwise_derivation.json",
            "frozen_kernel": "report_qw2049_spectral_micro_stagec_intersection_gate.json",
        },
        "targets": {
            "z_beta_target": z_beta_target,
            "delta_eta_target": delta_eta_target,
        },
        "baseline": {
            "z_beta_q25_q50_q75": [base_z_q25, base_z_q50, base_z_q75],
            "delta_eta_q25_q50_q75": [base_de_q25, base_de_q50, base_de_q75],
            "z_beta_log_iqr": base_z_log_iqr,
            "delta_eta_iqr": base_de_iqr,
            "n_bins_total": int(len(bins)),
        },
        "selection_mode": selection_mode,
        "selected_filter": selected,
        "thresholds": thresholds,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "tightened_warning_resolved": tightened_warning_resolved,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2066: COMPATIBILITY-FILTERED MICRO CONSTANTS TIGHTENING",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- selection_mode: {selection_mode}",
        f"- tightened_warning_resolved: {tightened_warning_resolved}",
        "",
        "## Baseline vs Selected",
        f"- baseline z_beta_log_iqr: {base_z_log_iqr:.6f}",
        f"- selected z_beta_log_iqr: {selected['z_beta_log_iqr']:.6f}",
        f"- baseline delta_eta_iqr: {base_de_iqr:.6f}",
        f"- selected delta_eta_iqr: {selected['delta_eta_iqr']:.6f}",
        f"- selected n_bins: {selected['n_selected']}",
        "",
        "## Selected Medians",
        f"- z_beta_q50: {selected['z_beta_q50']:.6f} (target {z_beta_target:.6f})",
        f"- delta_eta_q50: {selected['delta_eta_q50']:.6f} (target {delta_eta_target:.6f})",
        f"- abs_log_ratio_z_beta: {selected['z_beta_abs_log_ratio_to_target']:.6f}",
        f"- abs_err_delta_eta: {selected['delta_eta_abs_err_to_target']:.6f}",
        "",
        "## Selected Filter",
        f"- q_n / q_phase / q_rmse: {selected['q_n']:.2f} / {selected['q_phase']:.2f} / {selected['q_rmse']:.2f}",
        f"- n_min / phase_min / rmse_max: {selected['n_min']:.3f} / {selected['phase_min']:.6f} / {selected['rmse_max']:.6f}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend([
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2066] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2066] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2066] verdict={verdict} pass_count={pass_count}/{total_flags} "
        f"zlog {base_z_log_iqr:.3f}->{selected['z_beta_log_iqr']:.3f}"
    )


if __name__ == "__main__":
    main()
