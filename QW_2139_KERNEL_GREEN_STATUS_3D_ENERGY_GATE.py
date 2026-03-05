#!/usr/bin/env python3
"""
QW-2139: Kernel Green-status + 3D energy gate.

Purpose:
- close L14 to a rigorous PARTIAL state without overclaiming,
- separate three claims explicitly:
  1) classical local Green-function status,
  2) role of K as structural/index-space kernel in strict chain,
  3) 3D integrability class (absolute vs L2 / gradient energy).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2139_kernel_green_status_3d_energy_gate.json"
OUT_MD = ROOT / "RAPORT_QW2139_KERNEL_GREEN_STATUS_3D_ENERGY_GATE.md"


def kernel(r: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    den = 1.0 + beta * np.power(r, eta)
    return np.cos(omega * r + phi) / den


def dkernel(r: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    den = 1.0 + beta * np.power(r, eta)
    t = omega * r + phi
    term1 = -omega * np.sin(t) / den
    term2 = -(np.cos(t) * (beta * eta * np.power(r, eta - 1.0))) / (den * den)
    return term1 + term2


def l2_norm(x: np.ndarray, r: np.ndarray) -> float:
    return float(np.sqrt(np.trapz(x * x, r)))


def trapz_on_prefix(x: np.ndarray, r: np.ndarray, rmax: float) -> float:
    m = r <= rmax
    return float(np.trapz(x[m], r[m]))


def main() -> None:
    omega = 0.18575
    phi = 0.16250
    beta = 1.0
    eta = 1.8

    # Dense grid for stable radial diagnostics.
    r = np.linspace(1.0e-3, 200.0, 200_000)
    k = kernel(r, omega, phi, beta, eta)
    kp = dkernel(r, omega, phi, beta, eta)
    kpp = np.gradient(kp, r)

    # Radial Laplacian in 3D for r>0.
    lap_k = kpp + (2.0 / r) * kp

    # If K were harmonic away from source, lap_k should be ~0.
    norm_k = l2_norm(k, r)
    laplace_residual_ratio = l2_norm(lap_k, r) / max(norm_k, 1e-30)

    # Best-fit Helmholtz proxy: lap_k + lambda*K = 0.
    lam_best = -float(np.trapz(lap_k * k, r) / max(np.trapz(k * k, r), 1e-30))
    helm_res = lap_k + lam_best * k
    helmholtz_residual_ratio = l2_norm(helm_res, r) / max(norm_k, 1e-30)

    # 3D radial measures.
    abs_density = r * r * np.abs(k)
    l2_density = r * r * (k * k)
    grad_density = r * r * (kp * kp)

    rs: List[float] = [20.0, 40.0, 80.0, 120.0, 160.0, 200.0]
    abs_int = [trapz_on_prefix(abs_density, r, rv) for rv in rs]
    l2_int = [trapz_on_prefix(l2_density, r, rv) for rv in rs]
    grad_int = [trapz_on_prefix(grad_density, r, rv) for rv in rs]

    slope_abs = float(np.polyfit(np.log(np.array(rs)), np.log(np.array(abs_int)), 1)[0])
    expected_abs_slope = 3.0 - eta

    # Search for any dedicated constructive operator derivation artifact.
    # If absent, strict status cannot claim full Green-operator derivation.
    operator_artifacts = sorted(
        [p.name for p in ROOT.glob("report_qw*green*operator*.json")]
        + [p.name for p in ROOT.glob("report_qw*operator*green*.json")]
    )

    flags = {
        "frozen_kernel_release5_parameters_used": True,
        "laplace_green_condition_violated_away_from_source": bool(laplace_residual_ratio > 1e-1),
        "helmholtz_green_condition_violated_best_fit_lambda": bool(helmholtz_residual_ratio > 1e-1),
        "asymptotic_decay_not_classical_1_over_r": bool(abs(eta - 1.0) > 1e-6),
        "l1_3d_absolute_integral_diverges_for_eta_le_3": bool(eta <= 3.0),
        "l2_3d_energy_finite_for_eta_gt_1p5": bool(eta > 1.5),
        "gradient_energy_finite_for_eta_gt_0p5": bool(eta > 0.5),
        "numerical_abs_integral_growth_matches_3_minus_eta": bool(abs(slope_abs - expected_abs_slope) < 0.1),
        "no_constructive_green_operator_artifact_in_strict_chain": bool(len(operator_artifacts) == 0),
        "full_constructive_green_operator_derived_from_fin_action": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = len(flags)

    verdict = (
        "KERNEL_GREEN_STATUS_3D_ENERGY_GATE_PASS_PARTIAL_ROLE_CLARIFIED"
        if (
            flags["frozen_kernel_release5_parameters_used"]
            and flags["laplace_green_condition_violated_away_from_source"]
            and flags["helmholtz_green_condition_violated_best_fit_lambda"]
            and flags["asymptotic_decay_not_classical_1_over_r"]
            and flags["l1_3d_absolute_integral_diverges_for_eta_le_3"]
            and flags["l2_3d_energy_finite_for_eta_gt_1p5"]
            and flags["gradient_energy_finite_for_eta_gt_0p5"]
            and flags["numerical_abs_integral_growth_matches_3_minus_eta"]
            and flags["no_constructive_green_operator_artifact_in_strict_chain"]
        )
        else "KERNEL_GREEN_STATUS_3D_ENERGY_GATE_FAIL"
    )

    out: Dict = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_parameters": {
            "omega": omega,
            "phi": phi,
            "beta": beta,
            "eta": eta,
            "formula": "K(d)=cos(omega*d+phi)/(1+beta*d^eta)",
        },
        "green_status_checks": {
            "laplace_residual_ratio_l2": laplace_residual_ratio,
            "helmholtz_best_lambda": lam_best,
            "helmholtz_residual_ratio_l2": helmholtz_residual_ratio,
            "interpretation": (
                "K does not satisfy classical local Laplace/Helmholtz Green profile "
                "in tested radial diagnostics away from source."
            ),
        },
        "energy_3d_checks": {
            "radii_tested": rs,
            "integral_abs_r2_absK": abs_int,
            "integral_l2_r2_K2": l2_int,
            "integral_grad_r2_dK2": grad_int,
            "abs_growth_loglog_slope": slope_abs,
            "expected_abs_growth_slope_3_minus_eta": expected_abs_slope,
            "l2_tail_last_increment": float(l2_int[-1] - l2_int[-2]),
            "grad_tail_last_increment": float(grad_int[-1] - grad_int[-2]),
        },
        "operator_artifacts_found": operator_artifacts,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "classification": {
            "kernel_role_in_strict_chain": "STRUCTURAL_MIXING_KERNEL",
            "classical_local_green_function_claim": "NOT_SUPPORTED_IN_STRICT_CHAIN",
            "energy_statement_3d": (
                "ABSOLUTE_R2_WEIGHTED_INTEGRAL_DIVERGES_BUT_L2_AND_GRADIENT_ENERGY_ARE_FINITE_FOR_ETA_1P8"
            ),
        },
        "required_next_step": (
            "CONSTRUCT_EXPLICIT_OPERATOR_D_SUCH_THAT_DG_EQUALS_DELTA_OR_DECLARE_PERMANENT_NON_GREEN_ROLE_IN_CANONICAL_DOCS"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2139: KERNEL GREEN-STATUS + 3D ENERGY GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Kernel",
        f"- formula: `{out['kernel_parameters']['formula']}`",
        f"- omega={omega}, phi={phi}, beta={beta}, eta={eta}",
        "",
        "## Green-status diagnostics",
        f"- Laplace residual ratio: `{laplace_residual_ratio:.6f}`",
        f"- Helmholtz residual ratio (best lambda): `{helmholtz_residual_ratio:.6f}`",
        f"- best lambda: `{lam_best:.6f}`",
        "",
        "## 3D integrability diagnostics",
        f"- abs integral slope (log-log): `{slope_abs:.6f}` (expected `{expected_abs_slope:.6f}`)",
        f"- L2 tail last increment: `{out['energy_3d_checks']['l2_tail_last_increment']:.6e}`",
        f"- grad tail last increment: `{out['energy_3d_checks']['grad_tail_last_increment']:.6e}`",
        "",
        "## Scope boundary",
        "- Full constructive Green-operator derivation from FIN action: `False`",
        "",
    ]
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
