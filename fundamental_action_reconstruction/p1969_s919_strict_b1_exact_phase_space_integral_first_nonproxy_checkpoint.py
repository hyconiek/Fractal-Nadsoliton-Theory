#!/usr/bin/env python3
"""P1969 S919 strict B1 exact phase-space integral first non-proxy checkpoint.

This is a deliberately narrow continuation of P1864/P1863.  It replaces the
pointwise P1860 seed grid for the declared B1 projected Cutkosky polynomial by
an exact two-body massless phase-space integral, while preserving the no-false-
pass boundary: the integrated object is still the declared projected B1
polynomial kernel, not yet a full dressed graviton -> gauge_gauge amplitude.
"""

from __future__ import annotations

import json
import math
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy.integrate as integrate
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"

OUT = GEN / "p1969_s919_strict_b1_exact_phase_space_integral_first_nonproxy_checkpoint.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sympy_from_p1853(coeffs: dict[str, Any], key: str) -> sp.Expr:
    raw = (coeffs.get(key) or {}).get("symbolic")
    if raw is None:
        raise KeyError(f"missing symbolic coefficient {key}")
    return sp.sympify(raw, locals={"log": sp.log, "pi": sp.pi})


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")
    p1860 = load("p1860_s810_strict_b1_cutkosky_kernel_sample_and_pole_residue_table_checkpoint.json")
    p1863 = load("p1863_s813_strict_b1_projected_disc_uncertainty_certificate_checkpoint.json")

    coeffs = ((p1853.get("b1_symbolic_evaluation") or {}).get("evaluated_coefficients") or {})
    a_r2 = sympy_from_p1853(coeffs, "a_R2")
    a_ric2 = sympy_from_p1853(coeffs, "a_Ric2")

    # P1860 declared K_cut_sample(s, theta) = (a_R2 + a_Ric2*(1+theta^2))*s.
    # Here theta is promoted to x = cos(scattering angle) in [-1, 1], and the
    # azimuth is integrated explicitly with the massless two-body phase-space
    # convention dPi_2 = (1/(32*pi^2)) dphi dx.
    s, x, varphi = sp.symbols("s x varphi", positive=True, real=True)
    k_projected = s * (a_r2 + a_ric2 * (1 + x**2))
    d_pi_2_density = sp.Rational(1, 32) / sp.pi**2
    exact_integral = sp.simplify(
        sp.integrate(sp.integrate(d_pi_2_density * k_projected, (varphi, 0, 2 * sp.pi)), (x, -1, 1))
    )
    exact_expected = sp.simplify(s * (a_r2 / (8 * sp.pi) + a_ric2 / (6 * sp.pi)))
    exact_residual = sp.simplify(exact_integral - exact_expected)

    a_r2_f = float(sp.N(a_r2, 30))
    a_ric2_f = float(sp.N(a_ric2, 30))

    def integrand(phi_num: float, x_num: float, s_num: float) -> float:
        return (1.0 / (32.0 * math.pi**2)) * s_num * (a_r2_f + a_ric2_f * (1.0 + x_num * x_num))

    s_grid = np.array([0.5, 1.0, 2.0, 4.0, 8.0], dtype=float)
    rows: list[dict[str, Any]] = []
    for s_num in s_grid:
        numeric, error = integrate.dblquad(
            lambda phi_num, x_num: integrand(phi_num, x_num, float(s_num)),
            -1.0,
            1.0,
            lambda _x: 0.0,
            lambda _x: 2.0 * math.pi,
            epsabs=1e-14,
            epsrel=1e-14,
        )
        exact_num = float(sp.N(exact_integral.subs(s, float(s_num)), 30))
        abs_delta = abs(numeric - exact_num)
        # This is integration-error control only; it is not an amplitude-model
        # uncertainty for the still-missing dressed propagator/full vertex data.
        conservative_error = float(error + abs_delta + np.finfo(float).eps * max(1.0, abs(numeric)))
        lower_certified = numeric - conservative_error
        rows.append(
            {
                "s": float(s_num),
                "phase_space_integral_center": numeric,
                "sympy_exact_value": exact_num,
                "scipy_reported_abs_error": float(error),
                "abs_numeric_minus_exact": abs_delta,
                "conservative_integration_error": conservative_error,
                "lower_bound_from_integration_error": lower_certified,
                "positive_under_integration_error_only": lower_certified > 0.0,
            }
        )

    out = {
        "packet_id": "P1969",
        "stage_id": "S919",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "depends_on": {
            "p1853_present": "b1_symbolic_evaluation" in p1853,
            "p1860_present": "k_cut_sample_table" in p1860,
            "p1863_present": "projected_discontinuity_uncertainty_report" in p1863,
        },
        "kernel_scope": {
            "object_integrated": "P1860 declared projected B1 polynomial K_cut_sample promoted from point grid to angular variable x=cos(theta)",
            "formula": "K_proj(s,x)=s*(a_R2+a_Ric2*(1+x^2))",
            "phase_space_convention": "massless two-body dPi_2=(1/(32*pi^2))*dphi*dx, x in [-1,1], phi in [0,2*pi]",
            "strict_kernel_parameters_source": "P1853 evaluated coefficients from K_strict with omega=0.18575, phi=0.16250, beta=1.0, eta=1.8, alpha_geo_strict=4*ln(2)",
        },
        "symbolic_result": {
            "a_R2": str(a_r2),
            "a_Ric2": str(a_ric2),
            "integrand": str(k_projected),
            "phase_space_integral": str(exact_integral),
            "closed_form_check_target": str(exact_expected),
            "closed_form_residual": str(exact_residual),
            "closed_form_residual_is_zero": bool(exact_residual == 0),
        },
        "numeric_integral_table": {
            "s_grid": [float(v) for v in s_grid],
            "rows": rows,
            "all_positive_under_integration_error_only": all(r["positive_under_integration_error_only"] for r in rows),
            "minimum_lower_bound": min(r["lower_bound_from_integration_error"] for r in rows),
        },
        "strict_progress": {
            "replaces": "P1863 seed sigma=0.1*|disc_center| for this projected polynomial by exact symbolic integration plus scipy quadrature cross-check",
            "does_not_replace": [
                "full dressed graviton propagator pole computation",
                "full gauge-gauge vertex and polarization sum",
                "global RG-safe pole transport theorem",
                "background-family continuation beyond B1",
                "global unitarity or ToE closure",
            ],
        },
        "ur_link_theorem_candidate_boundary": {
            "positive_statement": "Within the P1860 projected B1 polynomial and declared phase-space convention, the integrated discontinuity is exactly s*(a_R2/(8*pi)+a_Ric2/(6*pi)) and is positive for the sampled physical s>0 corridor because the P1853 strict coefficients are positive.",
            "open_obstruction": "A genuine UR_link_theorem still needs dressed pole residues and RG flow covariance for the full amplitude, not just the projected polynomial backend.",
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": integrate.__version__ if hasattr(integrate, "__version__") else __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
        "false_pass_guard": "P1969 is an exact integral for the declared B1 projected polynomial only; it is not theorem-grade Cutkosky closure for the full strict theory.",
        "next_honest_step": "Compute the dressed graviton->gauge_gauge amplitude numerator and pole residues, then rerun this exact phase-space pipeline on the full polarization-summed integrand instead of the P1860 projected polynomial.",
        "lay_explanation": "Zamiast patrzeć na kilka punktów tabeli, policzyliśmy całkę po wszystkich kierunkach w najprostszym dwucząstkowym stanie końcowym. Wynik jest dodatni w tym ograniczonym modelu, ale pełny dowód wymaga jeszcze prawdziwej amplitudy i biegunów propagatora.",
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
