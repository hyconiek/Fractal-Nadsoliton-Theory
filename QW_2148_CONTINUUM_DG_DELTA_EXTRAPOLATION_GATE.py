#!/usr/bin/env python3
"""
QW-2148: Continuum DG=delta extrapolation gate (L14).

Purpose:
- quantify continuum-limit behavior from QW-2140/QW-2141 finite-grid evidence,
- keep strict distinction between numerical continuum support and action-level theorem,
- reduce ambiguity around the remaining L14 gap.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2148_continuum_dg_delta_extrapolation_gate.json"
OUT_MD = ROOT / "RAPORT_QW2148_CONTINUUM_DG_DELTA_EXTRAPOLATION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def monotone_nonincreasing(vals: List[float], tol: float = 1e-18) -> bool:
    return all(vals[i + 1] <= vals[i] + tol for i in range(len(vals) - 1))


def fit_error_limit(n_vals: List[float], e_vals: List[float], p: float) -> Tuple[float, float]:
    x = np.array([1.0 / (n**p) for n in n_vals], dtype=float)
    y = np.array(e_vals, dtype=float)
    # y = a*x + b, where b is extrapolated N->inf limit.
    a, b = np.polyfit(x, y, 1)
    yhat = a * x + b
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2)) if len(y) > 1 else 0.0
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(b), float(r2)


def main() -> None:
    r2140 = load("report_qw2140_kernel_inverse_finite_domain_gate.json")
    r2141 = load("report_qw2141_continuum_weak_distribution_proxy_gate.json")

    cases_2141 = sorted(r2141["cases"], key=lambda c: int(c["n_grid"]))
    n_grid = [int(c["n_grid"]) for c in cases_2141]
    e_reg = [float(c["max_abs_error_reg"]) for c in cases_2141]
    b_local = [float(c["max_boundary_sup_norm_local_only"]) for c in cases_2141]
    e_exact = [float(c["max_abs_error_exact"]) for c in cases_2141]

    reg_monotone = monotone_nonincreasing(e_reg)
    boundary_monotone = monotone_nonincreasing(b_local)

    b1, r2_1 = fit_error_limit(n_grid, e_reg, p=1.0)
    b2, r2_2 = fit_error_limit(n_grid, e_reg, p=2.0)
    if r2_2 >= r2_1:
        best_p = 2.0
        e_inf = b2
        best_r2 = r2_2
    else:
        best_p = 1.0
        e_inf = b1
        best_r2 = r2_1

    exact_inverse_errors = [
        float(c["exact_inverse_delta_metrics"]["delta_relative_l2_error"]) for c in r2140["cases"]
    ]
    max_exact_inverse_err = max(exact_inverse_errors)

    flags = {
        "finite_domain_inverse_constructive_available": bool(
            r2140["flags"]["constructive_finite_domain_inverse_operator_available"]
        ),
        "exact_inverse_delta_reconstruction_machine_precision": bool(max_exact_inverse_err <= 1e-12),
        "weak_distribution_proxy_regularized_error_monotone_with_volume": bool(reg_monotone),
        "boundary_aliasing_local_tests_monotone_down": bool(boundary_monotone),
        "extrapolated_continuum_error_limit_small": bool(e_inf <= 1e-6),
        "periodic_proxy_continuum_support_established": bool(
            reg_monotone and boundary_monotone and e_inf <= 1e-6
        ),
        "full_continuum_theorem_from_fin_action_completed": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = len(flags)

    verdict = (
        "CONTINUUM_DG_DELTA_EXTRAPOLATION_GATE_PASS_PARTIAL_ACTION_THEOREM_OPEN"
        if (
            flags["finite_domain_inverse_constructive_available"]
            and flags["exact_inverse_delta_reconstruction_machine_precision"]
            and flags["weak_distribution_proxy_regularized_error_monotone_with_volume"]
            and flags["boundary_aliasing_local_tests_monotone_down"]
            and flags["extrapolated_continuum_error_limit_small"]
            and flags["periodic_proxy_continuum_support_established"]
        )
        else "CONTINUUM_DG_DELTA_EXTRAPOLATION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2140": "report_qw2140_kernel_inverse_finite_domain_gate.json",
            "q2141": "report_qw2141_continuum_weak_distribution_proxy_gate.json",
        },
        "continuum_fit": {
            "n_grid": n_grid,
            "max_abs_error_reg": e_reg,
            "max_boundary_sup_norm_local_only": b_local,
            "max_abs_error_exact": e_exact,
            "best_power_p": best_p,
            "best_fit_r2": best_r2,
            "extrapolated_error_n_to_infinity": e_inf,
        },
        "finite_inverse_quality": {
            "max_exact_inverse_delta_relative_l2_error": max_exact_inverse_err,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "DERIVE_OPERATOR_D_AND_DISTRIBUTION_THEOREM_DK_EQUALS_DELTA_DIRECTLY_FROM_FIN_ACTION"
            if verdict.startswith("CONTINUUM_DG_DELTA_EXTRAPOLATION_GATE_PASS")
            else "REPAIR_CONVERGENCE_OR_PROXY_INPUTS_AND_RERUN_QW2148"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2148: CONTINUUM DG=DELTA EXTRAPOLATION GATE (L14)",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Continuum support summary",
        f"- n_grid: `{n_grid}`",
        f"- max_abs_error_reg: `{e_reg}`",
        f"- local boundary sup norms: `{b_local}`",
        f"- best fit power p: `{best_p}` (R2=`{best_r2:.6f}`)",
        f"- extrapolated error N->inf: `{e_inf:.6e}`",
        "",
        "## Boundary",
        "- Numerical continuum support is strengthened.",
        "- Full action-level theorem `DG=delta` remains open.",
        "",
    ]
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()

