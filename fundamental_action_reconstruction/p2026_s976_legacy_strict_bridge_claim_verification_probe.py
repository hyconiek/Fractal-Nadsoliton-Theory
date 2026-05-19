#!/usr/bin/env python3
"""P2026 S976: verify (or falsify) proposed legacy->strict bridge claim.

Guardrail-compliant probe:
- does NOT claim a bridge theorem,
- quantifies how far a proposed fractal-distance map can mimic strict kernel,
- exports explicit OPEN verdict with failure modes.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.optimize as so
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2026_s976_legacy_strict_bridge_claim_verification_probe.json"


def k_legacy(d: np.ndarray, alpha_geo: float, omega: float, phi: float, beta_tors: float) -> np.ndarray:
    return alpha_geo * np.cos(omega * d + phi) / (1.0 + beta_tors * d)


def k_strict(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d ** eta))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    d = np.linspace(1.0, 11.0, 200)

    # fixed strict lane reference used in existing p2025 seed
    omega_s, phi_s, beta_s, eta_s = 0.18575, 0.16250, 1.0, 1.8
    y_strict = k_strict(d, omega_s, phi_s, beta_s, eta_s)

    # legacy baseline (canonical amplitude present)
    alpha_geo = float(4.0 * np.log(2.0))

    # hypothesis: legacy can be warped into strict by d_eff=lambda*d^a and phase shift
    def objective(theta: np.ndarray) -> float:
        lam, a, dphi, beta_t = theta
        d_eff = lam * (d**a)
        y_leg = k_legacy(d_eff, alpha_geo, omega_s, phi_s + dphi, beta_t)
        # allow one scalar renormalization to avoid trivial amplitude mismatch
        c = float(np.dot(y_leg, y_strict) / max(np.dot(y_leg, y_leg), 1e-15))
        res = c * y_leg - y_strict
        return float(np.mean(res**2))

    bounds = [(0.1, 5.0), (0.01, 2.5), (-np.pi, np.pi), (1e-3, 10.0)]
    x0 = np.array([1.0, 1.0, 0.0, 1.0], dtype=float)
    res = so.minimize(objective, x0=x0, bounds=bounds, method="L-BFGS-B")
    lam, a, dphi, beta_t = [float(v) for v in res.x]
    d_eff = lam * (d**a)
    y_leg_best = k_legacy(d_eff, alpha_geo, omega_s, phi_s + dphi, beta_t)
    c_best = float(np.dot(y_leg_best, y_strict) / max(np.dot(y_leg_best, y_leg_best), 1e-15))
    y_fit = c_best * y_leg_best

    ss_res = float(np.sum((y_fit - y_strict) ** 2))
    ss_tot = float(np.sum((y_strict - np.mean(y_strict)) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1e-15)
    max_abs_gap = float(np.max(np.abs(y_fit - y_strict)))

    # symbolic non-identity witness: linear denominator cannot be identically equal to power eta=1.8
    d_sym, beta_t_sym, beta_s_sym = sp.symbols("d beta_t beta_s", positive=True, real=True)
    eta_sym = sp.Rational(9, 5)  # 1.8
    expr = sp.simplify((1 + beta_t_sym * d_sym) - (1 + beta_s_sym * d_sym**eta_sym))
    symbolic_identity_possible = bool(sp.simplify(expr) == 0)

    verdict = "OPEN_NONBRIDGE_EVIDENCE"
    if r2 > 0.99 and max_abs_gap < 1e-3 and symbolic_identity_possible:
        verdict = "CANDIDATE_BRIDGE_NEEDS_THEOREM"

    payload = {
        "schema_version": "p2026_s976_v1",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": verdict,
        "fit_summary": {
            "r2": float(r2),
            "max_abs_gap": max_abs_gap,
            "mse": float(ss_res / len(d)),
            "optimizer_success": bool(res.success),
            "optimizer_message": str(res.message),
        },
        "best_parameters": {
            "lambda": lam,
            "alpha_fractal": a,
            "delta_phi": dphi,
            "beta_tors_effective": beta_t,
            "amplitude_renorm_c": c_best,
        },
        "symbolic_check": {
            "linear_minus_fractional_denominator_expr": sp.sstr(expr),
            "identity_possible_without_new_axiom": symbolic_identity_possible,
        },
        "hard_limits": [
            "Numerical fit quality (even high R^2) is not a bridge theorem.",
            "No legacy-physical-role transfer to strict kernel is claimed.",
            "Selector closure (QW-2191) remains open.",
        ],
        "next_honest_step": "Build explicit theorem candidate that states necessary/sufficient assumptions under which d->lambda*d^alpha plus field normalization induces strict operational kernel, then test identifiability/uniqueness of those assumptions.",
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(str(OUT))


if __name__ == "__main__":
    main()

