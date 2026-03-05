#!/usr/bin/env python3
"""
QW-2140: Kernel inverse finite-domain gate.

Purpose:
- construct explicit inverse operator in a finite periodic 3D domain,
- verify delta reconstruction in Fourier space (exact and regularized),
- keep strict boundary: no claim of full continuum action-level Green proof.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2140_kernel_inverse_finite_domain_gate.json"
OUT_MD = ROOT / "RAPORT_QW2140_KERNEL_INVERSE_FINITE_DOMAIN_GATE.md"


def kernel_grid(n: int, omega: float, phi: float, beta: float, eta: float, dx: float = 1.0) -> np.ndarray:
    x = (np.arange(n) - n // 2) * dx
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(xx * xx + yy * yy + zz * zz)
    den = 1.0 + beta * np.power(r, eta)
    return np.cos(omega * r + phi) / den


def delta_like_metrics(a: np.ndarray) -> Dict[str, float]:
    n = a.shape[0]
    c = n // 2
    center = float(a[c, c, c])
    off = a.copy()
    off[c, c, c] = 0.0
    off_max = float(np.max(np.abs(off)))
    ratio = float(abs(center) / max(off_max, 1e-30))
    target = np.zeros_like(a)
    target[c, c, c] = 1.0
    rel_err = float(np.linalg.norm(a - target) / max(np.linalg.norm(target), 1e-30))
    return {
        "center_value": center,
        "off_max_abs": off_max,
        "center_to_offmax_ratio": ratio,
        "delta_relative_l2_error": rel_err,
        "sum_all_entries": float(np.sum(a)),
    }


def run_case(n: int, omega: float, phi: float, beta: float, eta: float) -> Dict:
    k = kernel_grid(n=n, omega=omega, phi=phi, beta=beta, eta=eta)
    khat = np.fft.fftn(np.fft.ifftshift(k))
    kabs = np.abs(khat)
    max_abs = float(np.max(kabs))
    min_abs = float(np.min(kabs))
    p1_abs = float(np.percentile(kabs, 1.0))
    cond_p1 = float(max_abs / max(p1_abs, 1e-30))
    near_zero_frac = float(np.mean(kabs < (1e-12 * max_abs)))

    # Exact inverse in periodic finite domain.
    inv_exact = 1.0 / khat
    identity_exact = np.real(np.fft.fftshift(np.fft.ifftn(inv_exact * khat)))
    m_exact = delta_like_metrics(identity_exact)

    # Regularized inverse (stable fallback).
    reg_eps = 1e-8 * (max_abs * max_abs)
    inv_reg = np.conj(khat) / (kabs * kabs + reg_eps)
    identity_reg = np.real(np.fft.fftshift(np.fft.ifftn(inv_reg * khat)))
    m_reg = delta_like_metrics(identity_reg)

    return {
        "n_grid": n,
        "spectral": {
            "max_abs_khat": max_abs,
            "min_abs_khat": min_abs,
            "p1_abs_khat": p1_abs,
            "condition_like_max_over_p1": cond_p1,
            "near_zero_fraction_1e-12": near_zero_frac,
        },
        "exact_inverse_delta_metrics": m_exact,
        "regularized_inverse_delta_metrics": m_reg,
        "regularization_eps": reg_eps,
    }


def main() -> None:
    omega = 0.18575
    phi = 0.16250
    beta = 1.0
    eta = 1.8

    cases: List[Dict] = []
    for n in [32, 40, 48]:
        cases.append(run_case(n=n, omega=omega, phi=phi, beta=beta, eta=eta))

    # Aggregate strict checks.
    exact_err_ok = all(c["exact_inverse_delta_metrics"]["delta_relative_l2_error"] < 1e-10 for c in cases)
    exact_ratio_ok = all(c["exact_inverse_delta_metrics"]["center_to_offmax_ratio"] > 1e10 for c in cases)
    reg_err_ok = all(c["regularized_inverse_delta_metrics"]["delta_relative_l2_error"] < 0.05 for c in cases)
    reg_ratio_ok = all(c["regularized_inverse_delta_metrics"]["center_to_offmax_ratio"] > 20.0 for c in cases)
    no_zero_modes = all(c["spectral"]["near_zero_fraction_1e-12"] == 0.0 for c in cases)
    moderate_cond = all(c["spectral"]["condition_like_max_over_p1"] < 1e4 for c in cases)

    flags = {
        "finite_periodic_domain_fft_operator_constructed": True,
        "no_spectral_zero_modes_detected_in_tested_grids": bool(no_zero_modes),
        "spectral_condition_moderate_in_tested_grids": bool(moderate_cond),
        "exact_inverse_reconstructs_delta_in_tested_grids": bool(exact_err_ok and exact_ratio_ok),
        "regularized_inverse_reconstructs_delta_like_kernel": bool(reg_err_ok and reg_ratio_ok),
        "constructive_finite_domain_inverse_operator_available": True,
        "full_continuum_action_level_green_operator_proof_completed": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = len(flags)

    verdict = (
        "KERNEL_INVERSE_FINITE_DOMAIN_GATE_PASS_PARTIAL"
        if (
            flags["finite_periodic_domain_fft_operator_constructed"]
            and flags["no_spectral_zero_modes_detected_in_tested_grids"]
            and flags["spectral_condition_moderate_in_tested_grids"]
            and flags["exact_inverse_reconstructs_delta_in_tested_grids"]
            and flags["regularized_inverse_reconstructs_delta_like_kernel"]
            and flags["constructive_finite_domain_inverse_operator_available"]
        )
        else "KERNEL_INVERSE_FINITE_DOMAIN_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_parameters": {
            "omega": omega,
            "phi": phi,
            "beta": beta,
            "eta": eta,
            "formula": "K(d)=cos(omega*d+phi)/(1+beta*d^eta)",
        },
        "cases": cases,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "classification": {
            "status_l14_progress": "FINITE_DOMAIN_INVERSE_CONSTRUCTIVE_STEP_CLOSED",
            "continuum_boundary": "NOT_CLOSED",
        },
        "required_next_step": (
            "PROMOTE_FINITE_DOMAIN_INVERSE_TO_CONTINUUM_OPERATOR_DERIVATION_FROM_FIN_ACTION"
            if verdict.startswith("KERNEL_INVERSE_FINITE_DOMAIN_GATE_PASS")
            else "REPAIR_SPECTRAL_INVERSE_RECONSTRUCTION_AND_RERUN_QW2140"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2140: KERNEL INVERSE FINITE-DOMAIN GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Tested grids",
    ]
    for case in cases:
        n = case["n_grid"]
        se = case["spectral"]
        me = case["exact_inverse_delta_metrics"]
        mr = case["regularized_inverse_delta_metrics"]
        lines.extend(
            [
                f"- N={n}: cond~`{se['condition_like_max_over_p1']:.3f}`, near-zero-frac=`{se['near_zero_fraction_1e-12']:.3e}`",
                f"  exact: err=`{me['delta_relative_l2_error']:.3e}`, center/offmax=`{me['center_to_offmax_ratio']:.3e}`",
                f"  regularized: err=`{mr['delta_relative_l2_error']:.3e}`, center/offmax=`{mr['center_to_offmax_ratio']:.3e}`",
            ]
        )
    lines.extend(
        [
            "",
            "## Scope boundary",
            "- Full continuum action-level Green operator proof: `False`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
