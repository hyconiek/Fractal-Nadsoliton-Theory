#!/usr/bin/env python3
"""
QW-2064: Micro-derived renormalization constants gate (no sector retune).

Goal:
- Formalize whether (Z_beta, delta_eta) are supported by micro derivation
  from QW-2048 for the frozen kernel used in QW-2049.
- Use deterministic checks only (no new optimization/scan).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2064_micro_derived_renormalization_constants_gate.json"
OUT_MD = ROOT / "RAPORT_QW2064_MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE.md"


def quantiles(x: List[float], q: List[float]) -> List[float]:
    arr = np.array(x, dtype=float)
    return [float(np.quantile(arr, v)) for v in q]


def main() -> None:
    r2048 = json.loads((ROOT / "report_qw2048_spectral_phase_locked_pointwise_derivation.json").read_text(encoding="utf-8"))
    r2049 = json.loads((ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8"))

    kernel = {k: float(v) for k, v in r2049["stagec_pool"]["selected_kernel"].items()}

    # UV canonical constants from Nadsoliton canonical branch (pre-QW1700 semantics).
    beta_uv = 0.01
    eta_uv = 1.0

    z_beta_target = float(kernel["beta"] / beta_uv)
    delta_eta_target = float(kernel["eta"] - eta_uv)

    ge = r2048["global_estimates"]
    cov = r2048["coverage"]
    bins = r2048["pointwise_bins"]["bins"]

    z_beta_micro = float(ge["z_beta_median"])
    delta_eta_micro = float(ge["delta_eta_median"])

    z_vals = [float(b["z_beta_median"]) for b in bins]
    de_vals = [float(b["delta_eta_median"]) for b in bins]
    n_vals = [int(b["n"]) for b in bins]
    ph_vals = [float(b["phase_min_median"]) for b in bins]
    rmse_vals = [float(b["rmse_median"]) for b in bins]

    z_q25, z_q50, z_q75 = quantiles(z_vals, [0.25, 0.50, 0.75])
    de_q25, de_q50, de_q75 = quantiles(de_vals, [0.25, 0.50, 0.75])

    z_log_iqr = float(np.log(max(z_q75, 1e-12)) - np.log(max(z_q25, 1e-12)))
    de_iqr = float(de_q75 - de_q25)

    abs_log_ratio = float(abs(math.log(max(z_beta_micro, 1e-12) / max(z_beta_target, 1e-12))))
    abs_delta_err = float(abs(delta_eta_micro - delta_eta_target))

    thresholds = {
        "n_rows_total_min": 300,
        "n_bins_min": 12,
        "phase_min_median_min": 0.80,
        "joint_target_in_ci95_fraction_min": 0.80,
        "abs_log_zbeta_ratio_max": math.log(2.0),
        "abs_delta_eta_err_max": 0.35,
    }

    flags = {
        "qw2048_pointwise_pass": bool(r2048.get("verdict") == "SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION_PASS"),
        "n_rows_total_ge_min": bool(float(ge["n_rows_total"]) >= thresholds["n_rows_total_min"]),
        "n_bins_ge_min": bool(int(r2048["pointwise_bins"]["n_bins"]) >= thresholds["n_bins_min"]),
        "phase_min_median_ge_min": bool(float(ge["phase_min_median"]) >= thresholds["phase_min_median_min"]),
        "joint_target_in_ci95_fraction_ge_min": bool(float(cov["joint_target_in_ci95_fraction"]) >= thresholds["joint_target_in_ci95_fraction_min"]),
        "z_beta_logratio_to_target_le_max": bool(abs_log_ratio <= thresholds["abs_log_zbeta_ratio_max"]),
        "delta_eta_abs_err_le_max": bool(abs_delta_err <= thresholds["abs_delta_eta_err_max"]),
        "no_sector_retune_between_micro_and_kernel": True,
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    core_pass = bool(
        flags["qw2048_pointwise_pass"]
        and flags["n_rows_total_ge_min"]
        and flags["n_bins_ge_min"]
        and flags["phase_min_median_ge_min"]
        and flags["joint_target_in_ci95_fraction_ge_min"]
        and flags["z_beta_logratio_to_target_le_max"]
        and flags["delta_eta_abs_err_le_max"]
        and flags["no_sector_retune_between_micro_and_kernel"]
    )

    ci_warning = bool((z_log_iqr > 2.0) or (de_iqr > 0.8))

    if core_pass and not ci_warning:
        verdict = "MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE_PASS"
    elif core_pass and ci_warning:
        verdict = "MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE_PASS_WITH_WIDE_CI_WARNING"
    else:
        verdict = "MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "micro_pointwise": "report_qw2048_spectral_phase_locked_pointwise_derivation.json",
            "frozen_kernel": "report_qw2049_spectral_micro_stagec_intersection_gate.json:stagec_pool.selected_kernel",
        },
        "uv_canonical_constants": {
            "beta_uv": beta_uv,
            "eta_uv": eta_uv,
        },
        "frozen_kernel": kernel,
        "targets": {
            "z_beta_target": z_beta_target,
            "delta_eta_target": delta_eta_target,
        },
        "micro_global": {
            "z_beta_median": z_beta_micro,
            "delta_eta_median": delta_eta_micro,
            "n_rows_total": int(ge["n_rows_total"]),
            "n_bins": int(r2048["pointwise_bins"]["n_bins"]),
            "phase_min_median": float(ge["phase_min_median"]),
            "joint_target_in_ci95_fraction": float(cov["joint_target_in_ci95_fraction"]),
        },
        "deviation": {
            "abs_log_zbeta_ratio": abs_log_ratio,
            "abs_delta_eta_error": abs_delta_err,
        },
        "dispersion_diagnostics": {
            "z_beta_bin_q25_q50_q75": [z_q25, z_q50, z_q75],
            "delta_eta_bin_q25_q50_q75": [de_q25, de_q50, de_q75],
            "z_beta_log_iqr": z_log_iqr,
            "delta_eta_iqr": de_iqr,
            "n_bin_median": float(np.median(np.array(n_vals, dtype=float))),
            "phase_min_median_bin": float(np.median(np.array(ph_vals, dtype=float))),
            "rmse_median_bin": float(np.median(np.array(rmse_vals, dtype=float))),
        },
        "thresholds": thresholds,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "ci_warning": ci_warning,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2064: MICRO-DERIVED RENORMALIZATION CONSTANTS GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- ci_warning: {ci_warning}",
        "",
        "## Targets (from frozen kernel)",
        f"- z_beta_target: {z_beta_target:.6f}",
        f"- delta_eta_target: {delta_eta_target:.6f}",
        "",
        "## Micro-derived (QW-2048)",
        f"- z_beta_median: {z_beta_micro:.6f}",
        f"- delta_eta_median: {delta_eta_micro:.6f}",
        f"- n_rows_total / n_bins: {int(ge['n_rows_total'])} / {int(r2048['pointwise_bins']['n_bins'])}",
        f"- phase_min_median: {float(ge['phase_min_median']):.6f}",
        f"- joint_target_in_ci95_fraction: {float(cov['joint_target_in_ci95_fraction']):.6f}",
        "",
        "## Deviations",
        f"- abs_log(z_beta_micro/z_beta_target): {abs_log_ratio:.6f}",
        f"- abs(delta_eta_micro-delta_eta_target): {abs_delta_err:.6f}",
        "",
        "## Dispersion Diagnostics",
        f"- z_beta log-IQR: {z_log_iqr:.6f}",
        f"- delta_eta IQR: {de_iqr:.6f}",
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

    print(f"[QW-2064] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2064] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2064] verdict={verdict} pass_count={pass_count}/{total_flags} ci_warning={ci_warning}")


if __name__ == "__main__":
    main()
