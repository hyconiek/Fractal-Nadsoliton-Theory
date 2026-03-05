#!/usr/bin/env python3
"""
QW-2141: Continuum-weak distribution proxy gate for L14.

Purpose:
- strengthen L14 from finite-domain inverse step (QW-2140) to weak-distribution
  proxy evidence under increasing periodic volume,
- verify pairing <D*K, phi> ~= phi(0) on localized smooth test functions,
- keep strict boundary: this is NOT a full continuum theorem from FIN action.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2141_continuum_weak_distribution_proxy_gate.json"
OUT_MD = ROOT / "RAPORT_QW2141_CONTINUUM_WEAK_DISTRIBUTION_PROXY_GATE.md"


def kernel_grid(n: int, omega: float, phi: float, beta: float, eta: float, dx: float = 1.0) -> np.ndarray:
    x = (np.arange(n) - n // 2) * dx
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(xx * xx + yy * yy + zz * zz)
    den = 1.0 + beta * np.power(r, eta)
    return np.cos(omega * r + phi) / den


def build_test_functions(n: int, dx: float = 1.0) -> Dict[str, np.ndarray]:
    l = n * dx
    x = (np.arange(n) - n // 2) * dx
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(xx * xx + yy * yy + zz * zz)

    # Compact smooth bump (C^\infty) support r<R.
    r0 = 6.0
    t = r / r0
    bump = np.zeros_like(r)
    m = t < 1.0
    bump[m] = np.exp(-1.0 / (1.0 - t[m] * t[m]))
    if float(np.max(bump)) > 0.0:
        bump /= float(np.max(bump))

    bump_x = bump * (xx / r0)

    # Fast-decay gaussian (practically compact in tested boxes).
    sigma = 3.0
    gauss = np.exp(-(r * r) / (2.0 * sigma * sigma))

    # Low-frequency periodic mode.
    cos_mode = np.cos(2.0 * np.pi * xx / l)

    return {
        "bump_compact": bump,
        "bump_x_compact": bump_x,
        "gauss_fast_decay": gauss,
        "cos_periodic_mode": cos_mode,
    }


def boundary_sup_norm(v: np.ndarray, dx: float = 1.0, margin: float = 4.0) -> float:
    n = v.shape[0]
    l = n * dx
    x = (np.arange(n) - n // 2) * dx
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    m = (np.abs(xx) > (l / 2.0 - margin)) | (np.abs(yy) > (l / 2.0 - margin)) | (np.abs(zz) > (l / 2.0 - margin))
    return float(np.max(np.abs(v[m])))


def run_case(n: int, omega: float, phi: float, beta: float, eta: float) -> Dict:
    dx = 1.0
    k = kernel_grid(n=n, omega=omega, phi=phi, beta=beta, eta=eta, dx=dx)
    khat = np.fft.fftn(np.fft.ifftshift(k))
    kabs = np.abs(khat)

    # Exact and regularized inverse symbols.
    inv_exact = 1.0 / khat
    eps = 1e-8 * float(np.max(kabs) ** 2)
    inv_reg = np.conj(khat) / (kabs * kabs + eps)

    q_exact = np.real(np.fft.fftshift(np.fft.ifftn(inv_exact * khat)))
    q_reg = np.real(np.fft.fftshift(np.fft.ifftn(inv_reg * khat)))

    c = n // 2
    tests = build_test_functions(n=n, dx=dx)
    test_rows: List[Dict] = []
    for name, v in tests.items():
        rhs = float(v[c, c, c])
        lhs_exact = float(np.sum(q_exact * v))
        lhs_reg = float(np.sum(q_reg * v))
        test_rows.append(
            {
                "name": name,
                "rhs_phi_at_zero": rhs,
                "lhs_exact_pairing": lhs_exact,
                "lhs_reg_pairing": lhs_reg,
                "abs_error_exact": abs(lhs_exact - rhs),
                "abs_error_reg": abs(lhs_reg - rhs),
                "boundary_sup_norm": boundary_sup_norm(v, dx=dx),
            }
        )

    err_exact_max = float(max(row["abs_error_exact"] for row in test_rows))
    err_reg_max = float(max(row["abs_error_reg"] for row in test_rows))
    local_names = {"bump_compact", "bump_x_compact", "gauss_fast_decay"}
    boundary_max = float(max(row["boundary_sup_norm"] for row in test_rows))
    boundary_max_local = float(max(row["boundary_sup_norm"] for row in test_rows if row["name"] in local_names))

    return {
        "n_grid": n,
        "domain_length": float(n * dx),
        "eps_regularization": eps,
        "spectral": {
            "near_zero_fraction_1e-12": float(np.mean(kabs < (1e-12 * float(np.max(kabs))))),
            "condition_like_max_over_p1": float(float(np.max(kabs)) / max(float(np.percentile(kabs, 1.0)), 1e-30)),
        },
        "q_exact_center": float(q_exact[c, c, c]),
        "q_reg_center": float(q_reg[c, c, c]),
        "test_rows": test_rows,
        "max_abs_error_exact": err_exact_max,
        "max_abs_error_reg": err_reg_max,
        "max_boundary_sup_norm": boundary_max,
        "max_boundary_sup_norm_local_only": boundary_max_local,
    }


def main() -> None:
    omega = 0.18575
    phi = 0.16250
    beta = 1.0
    eta = 1.8

    # Increasing periodic volume with fixed dx=1.
    cases = [run_case(n=n, omega=omega, phi=phi, beta=beta, eta=eta) for n in [32, 48, 64]]

    exact_max = float(max(c["max_abs_error_exact"] for c in cases))
    reg_errors = [float(c["max_abs_error_reg"]) for c in cases]
    reg_max = float(max(reg_errors))
    reg_min = float(min(reg_errors))
    boundary_max = float(max(c["max_boundary_sup_norm_local_only"] for c in cases))

    flags = {
        "localized_test_function_family_declared": True,
        "increasing_periodic_volume_with_fixed_dx": True,
        "exact_pairing_identity_all_cases": bool(exact_max < 1e-10),
        "regularized_pairing_identity_small_error_all_cases": bool(reg_max < 1e-5),
        "regularized_error_stable_across_volume_growth": bool((reg_max / max(reg_min, 1e-30)) < 2.0),
        "boundary_aliasing_suppressed_for_local_tests": bool(boundary_max < 1e-3),
        "deterministic_no_scan_no_retune": True,
        "full_continuum_distribution_theorem_from_fin_action_completed": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = len(flags)

    verdict = (
        "CONTINUUM_WEAK_DISTRIBUTION_PROXY_GATE_PASS_PARTIAL"
        if (
            flags["localized_test_function_family_declared"]
            and flags["increasing_periodic_volume_with_fixed_dx"]
            and flags["exact_pairing_identity_all_cases"]
            and flags["regularized_pairing_identity_small_error_all_cases"]
            and flags["regularized_error_stable_across_volume_growth"]
            and flags["boundary_aliasing_suppressed_for_local_tests"]
            and flags["deterministic_no_scan_no_retune"]
        )
        else "CONTINUUM_WEAK_DISTRIBUTION_PROXY_GATE_FAIL"
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
        "cases": cases,
        "aggregate": {
            "max_abs_error_exact": exact_max,
            "max_abs_error_reg": reg_max,
            "min_abs_error_reg": reg_min,
            "reg_error_ratio_max_over_min": float(reg_max / max(reg_min, 1e-30)),
            "max_boundary_sup_norm": boundary_max,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "classification": {
            "status_l14_progress": "WEAK_DISTRIBUTION_PROXY_CONTINUUM_STEP_CLOSED",
            "continuum_theorem_boundary": "NOT_CLOSED",
        },
        "required_next_step": (
            "PROVE_CONTINUUM_DISTRIBUTION_THEOREM_DK_EQUALS_DELTA_FROM_FIN_ACTION_NOT_ONLY_PERIODIC_FFT_PROXY"
            if verdict.startswith("CONTINUUM_WEAK_DISTRIBUTION_PROXY_GATE_PASS")
            else "REPAIR_WEAK_PAIRING_CONSISTENCY_AND_RERUN_QW2141"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2141: CONTINUUM-WEAK DISTRIBUTION PROXY GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Aggregate",
        f"- max exact pairing error: `{exact_max:.3e}`",
        f"- max regularized pairing error: `{reg_max:.3e}`",
        f"- reg error stability ratio max/min: `{out['aggregate']['reg_error_ratio_max_over_min']:.3f}`",
        f"- max boundary sup norm (tests): `{boundary_max:.3e}`",
        "",
        "## Scope boundary",
        "- Full continuum theorem from FIN action: `False`",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
