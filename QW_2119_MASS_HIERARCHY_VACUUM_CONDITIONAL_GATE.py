#!/usr/bin/env python3
"""
QW-2119: Mass-hierarchy exponential test + conditional vacuum closure gate.

Purpose:
- test whether current derived mass rows follow exponential hierarchy structure,
- quantify vacuum-closure requirement implied by lambda_min(K_total),
- keep strict separation: hierarchy evidence vs missing absolute scalar scale input.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2119_mass_hierarchy_vacuum_conditional_gate.json"
OUT_MD = ROOT / "RAPORT_QW2119_MASS_HIERARCHY_VACUUM_CONDITIONAL_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def kernel_value(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return float(np.cos(omega * d + phi) / (1.0 + beta * (d**eta)))


def build_ktotal_matrix(n: int, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    m = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                m[i, j] = 0.0
                continue
            d = min(abs(i - j), n - abs(i - j))
            m[i, j] = kernel_value(float(d), omega, phi, beta, eta)
    return 0.5 * (m + m.T)


def linear_log_fit(x: np.ndarray, y_mass: np.ndarray) -> Tuple[float, float, float]:
    yy = np.log(np.clip(y_mass, 1e-300, None))
    a = np.vstack([np.ones_like(x), x]).T
    coef, _, _, _ = np.linalg.lstsq(a, yy, rcond=None)
    yhat = a @ coef
    ss_res = float(np.sum((yy - yhat) ** 2))
    ss_tot = float(np.sum((yy - np.mean(yy)) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1e-300)
    intercept = float(coef[0])
    slope = float(coef[1])
    return intercept, slope, float(r2)


def main() -> None:
    r2063 = load_json("report_qw2063_derivational_reconstruction_shared_flavor_basis.json")
    r2049 = load_json("report_qw2049_spectral_micro_stagec_intersection_gate.json")
    r2118 = load_json("report_qw2118_ktotal_spectral_tripartition_gate.json")

    rows = r2063["metrics"]["mass"]["rows"]
    q_eff = np.array([float(r["q_eff"]) for r in rows], dtype=float)
    m_pred = np.array([float(r["pred_mev"]) for r in rows], dtype=float)
    m_exp = np.array([float(r["exp_mev"]) for r in rows], dtype=float)

    pred_intercept, pred_slope, pred_r2 = linear_log_fit(q_eff, m_pred)
    exp_intercept, exp_slope, exp_r2 = linear_log_fit(q_eff, m_exp)
    slope_rel_diff_pct = float(abs(pred_slope - exp_slope) / max(abs(exp_slope), 1e-15) * 100.0)

    ratio_pred = float(np.max(m_pred) / max(np.min(m_pred), 1e-300))
    ratio_exp = float(np.max(m_exp) / max(np.min(m_exp), 1e-300))
    ratio_rel_err_pct = float(abs(ratio_pred - ratio_exp) / max(abs(ratio_exp), 1e-15) * 100.0)

    kernel = {k: float(v) for k, v in r2049["stagec_pool"]["selected_kernel"].items()}
    m_ktotal = build_ktotal_matrix(
        12, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"]
    )
    lam_min = float(np.min(np.linalg.eigvalsh(m_ktotal)))
    required_shift = float(max(0.0, -lam_min + 1e-9))

    # Crucial strict distinction:
    # we can quantify required shift, but cannot close absolute vacuum stability
    # without explicit scalar-sector scale (vev / mass-floor calibration in same units).
    has_absolute_scalar_scale_input = False
    vacuum_closed_without_scale_input = False

    flags = {
        "hierarchy_pred_r2_ge_0p999": bool(pred_r2 >= 0.999),
        "hierarchy_exp_r2_ge_0p98": bool(exp_r2 >= 0.98),
        "hierarchy_slope_rel_diff_le_10pct": bool(slope_rel_diff_pct <= 10.0),
        "hierarchy_maxmin_ratio_rel_err_le_30pct": bool(ratio_rel_err_pct <= 30.0),
        "required_mass_shift_quantified": bool(required_shift >= 0.0),
        "vacuum_closed_without_scale_input": bool(vacuum_closed_without_scale_input),
        "tripartition_gate_consistent": bool(
            str(r2118.get("verdict", "")).startswith("KTOTAL_SPECTRAL_TRIPARTITION_GATE_PASS")
        ),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    hierarchy_pass = bool(
        flags["hierarchy_pred_r2_ge_0p999"]
        and flags["hierarchy_exp_r2_ge_0p98"]
        and flags["hierarchy_slope_rel_diff_le_10pct"]
        and flags["hierarchy_maxmin_ratio_rel_err_le_30pct"]
        and flags["required_mass_shift_quantified"]
    )

    verdict = (
        "MASS_HIERARCHY_PASS_VACUUM_CLOSURE_CONDITIONAL_ON_SCALE_INPUT"
        if hierarchy_pass and (not vacuum_closed_without_scale_input)
        else "MASS_HIERARCHY_AND_VACUUM_GATE_FAIL_OR_INCOMPLETE"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "mass_rows": "report_qw2063_derivational_reconstruction_shared_flavor_basis.json:metrics.mass.rows",
            "kernel": "report_qw2049_spectral_micro_stagec_intersection_gate.json:stagec_pool.selected_kernel",
            "tripartition_gate": "report_qw2118_ktotal_spectral_tripartition_gate.json",
        },
        "hierarchy_fit": {
            "pred": {
                "log_fit_intercept": pred_intercept,
                "log_fit_slope": pred_slope,
                "r2": pred_r2,
            },
            "exp": {
                "log_fit_intercept": exp_intercept,
                "log_fit_slope": exp_slope,
                "r2": exp_r2,
            },
            "slope_rel_diff_pct": slope_rel_diff_pct,
            "maxmin_ratio_pred": ratio_pred,
            "maxmin_ratio_exp": ratio_exp,
            "maxmin_ratio_rel_err_pct": ratio_rel_err_pct,
        },
        "vacuum_condition": {
            "lambda_min_ktotal": lam_min,
            "required_uniform_mass_shift_ge": required_shift,
            "has_absolute_scalar_scale_input": has_absolute_scalar_scale_input,
            "vacuum_closed_without_scale_input": vacuum_closed_without_scale_input,
            "note": (
                "Absolute vacuum closure requires explicit scalar-sector scale in units compatible with K_total."
            ),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROVIDE_EXPLICIT_SCALAR_SCALE_AND_CLOSE_VACUUM_STABILITY_WITH_QW2120"
            if verdict == "MASS_HIERARCHY_PASS_VACUUM_CLOSURE_CONDITIONAL_ON_SCALE_INPUT"
            else "REVIEW_HIERARCHY_OR_VACUUM_ASSUMPTIONS_AND_RERUN_QW2119"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2119: MASS HIERARCHY + VACUUM CONDITIONAL GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Exponential hierarchy fit",
        f"- pred r2 / slope: `{pred_r2:.6f}` / `{pred_slope:.9f}`",
        f"- exp  r2 / slope: `{exp_r2:.6f}` / `{exp_slope:.9f}`",
        f"- slope rel diff pct: `{slope_rel_diff_pct:.6f}`",
        f"- max/min ratio pred vs exp rel err pct: `{ratio_rel_err_pct:.6f}`",
        "",
        "## Vacuum closure condition",
        f"- lambda_min(K_total): `{lam_min:.9f}`",
        f"- required uniform mass shift >= `{required_shift:.9f}`",
        f"- absolute scalar scale input available: `{has_absolute_scalar_scale_input}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2119] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2119] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2119] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

