#!/usr/bin/env python3
"""P1981 S931 strict ADM/Bianchi-I R^2 lapse-variation obligation.

This extends P1980 from the Einstein-Hilbert bulk term to the first
curvature-squared term, R^2.  It keeps the lapse N(t), includes the Ndot
dependence of the Bianchi-I Ricci scalar, and computes the Euler-Lagrange lapse
operator for the anisotropic R^2 density difference.

Scope: R^2 minisuperspace lapse equation only.  It does not cover Ricci^2,
Riemann^2, Gauss-Bonnet, spatial equations, or global background independence.
"""

from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"
OUT = GEN / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1980 = load("p1980_s930_strict_adm_bianchi_eh_lapse_variation_witness.json")

    N, Nd, Ndd = sp.symbols("N Nd Ndd", positive=True, real=True)
    V, H, Hd, Hdd = sp.symbols("V H Hd Hdd", real=True)
    Q, Qd = sp.symbols("Q Qd", real=True)
    sigma1, sigma2 = sp.symbols("sigma1 sigma2", real=True)

    q_shear = sp.factor(sigma1**2 + sigma1 * sigma2 + sigma2**2)
    q_matrix = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    scipy_eigs = la.eigvalsh(np.array(q_matrix.tolist(), dtype=float))
    q_positive_definite = all(float(ev) > 0 for ev in scipy_eigs)

    # Ricci scalar for diagonal Bianchi-I with lapse N(t):
    # R = 2/N^2 [sum_i dot(h_i) + sum_i h_i^2 + sum_{i<j} h_i h_j
    #            - (Ndot/N) sum_i h_i], with h_i = dot(a_i)/a_i.
    # Under trace-free shear, sum_i h_i = 3H and the anisotropic excess is Q.
    A_frw = sp.factor(3 * Hd + 6 * H**2 - 3 * H * Nd / N)
    ricci_frw = sp.factor(2 * A_frw / N**2)
    ricci_bianchi = sp.factor(2 * (A_frw + Q) / N**2)
    ricci_difference = sp.factor(ricci_bianchi - ricci_frw)

    # R^2 density with the volume V=a1*a2*a3.  The lapse Euler-Lagrange operator
    # must include Ndot because R contains Ndot/N.
    r2_density_difference = sp.factor(N * V * (ricci_bianchi**2 - ricci_frw**2))

    total_derivative_rules = {
        N: Nd,
        Nd: Ndd,
        V: 3 * H * V,
        H: Hd,
        Hd: Hdd,
        Q: Qd,
    }

    def total_dt(expr: sp.Expr) -> sp.Expr:
        return sp.factor(sum(sp.diff(expr, var) * dvar for var, dvar in total_derivative_rules.items()))

    d_l_d_n = sp.factor(sp.diff(r2_density_difference, N))
    d_l_d_nd = sp.factor(sp.diff(r2_density_difference, Nd))
    dt_d_l_d_nd = sp.factor(total_dt(d_l_d_nd))
    lapse_el_difference = sp.factor(d_l_d_n - dt_d_l_d_nd)
    expected_lapse_el_difference = sp.factor(-12 * V * (6 * H**2 * Q - 2 * H * Qd + 4 * Hd * Q + Q**2) / N**4)
    el_matches_expected = sp.simplify(lapse_el_difference - expected_lapse_el_difference) == 0
    isotropic_limit_zero = sp.simplify(lapse_el_difference.subs({Q: 0, Qd: 0})) == 0
    contains_qd = bool(lapse_el_difference.has(Qd))
    contains_q_squared = bool(lapse_el_difference.has(Q**2))

    sample_points = [
        {N: sp.Rational(1, 1), V: sp.Rational(1, 1), H: sp.Rational(1, 1), Hd: sp.Rational(1, 10), Q: sp.Rational(3, 400), Qd: sp.Rational(1, 500)},
        {N: sp.Rational(5, 4), V: sp.Rational(7, 5), H: sp.Rational(3, 2), Hd: sp.Rational(-1, 20), Q: sp.Rational(7, 100), Qd: sp.Rational(-3, 1000)},
        {N: sp.Rational(4, 3), V: sp.Rational(11, 10), H: sp.Rational(4, 5), Hd: sp.Rational(1, 25), Q: sp.Rational(3, 256), Qd: sp.Rational(1, 640)},
        {N: sp.Rational(9, 10), V: sp.Rational(6, 5), H: sp.Rational(5, 4), Hd: sp.Rational(-1, 30), Q: sp.Rational(247, 3969), Qd: sp.Rational(-2, 1701)},
    ]
    rows: list[dict[str, Any]] = []
    for idx, point in enumerate(sample_points):
        expr_val = sp.simplify(lapse_el_difference.subs(point))
        expected_val = sp.simplify(expected_lapse_el_difference.subs(point))
        rows.append(
            {
                "sample_id": f"r2_lapse_variation_sample_{idx}",
                "input": {str(k): str(v) for k, v in point.items()},
                "lapse_EL_difference": str(expr_val),
                "expected_lapse_EL_difference": str(expected_val),
                "absolute_float_value": float(abs(sp.N(expr_val, 30))),
                "matches_expected": bool(sp.simplify(expr_val - expected_val) == 0),
            }
        )

    p1980_pass = p1980.get("result_kind") == "PASS_EH_LAPSE_SHEAR_TERM_DERIVED"
    all_samples_match = all(row["matches_expected"] for row in rows)
    theorem_pass = bool(p1980_pass and el_matches_expected and isotropic_limit_zero and contains_qd and contains_q_squared and q_positive_definite and all_samples_match)

    out = {
        "ledger_id": "P1981_S931_STRICT_ADM_BIANCHI_R2_LAPSE_VARIATION_OBLIGATION",
        "packet_id": "P1981",
        "stage_id": "S931",
        "produced_by": "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1980_eh_lapse_witness_present": p1980_pass,
        },
        "ricci_scalar_setup": {
            "coordinate_rates": "h_i = dot(a_i)/a_i, sum_i sigma_i=0",
            "A_FRW": str(A_frw),
            "R_FRW": str(ricci_frw),
            "R_BianchiI_tracefree_shear": str(ricci_bianchi),
            "R_BianchiI_minus_R_FRW": str(ricci_difference),
            "Q_shear_symbolic": str(q_shear),
            "Q_matrix": str(q_matrix),
            "Q_eigenvalues_exact": ["1/2", "3/2"],
            "Q_eigenvalues_scipy": [float(ev) for ev in scipy_eigs],
            "Q_positive_definite": q_positive_definite,
        },
        "r2_lapse_euler_operator": {
            "density_difference_NV_R2": str(r2_density_difference),
            "dL_dN": str(d_l_d_n),
            "dL_dNdot": str(d_l_d_nd),
            "Dt_dL_dNdot": str(dt_d_l_d_nd),
            "EL_N_difference": str(lapse_el_difference),
            "expected_EL_N_difference": str(expected_lapse_el_difference),
            "EL_matches_expected": el_matches_expected,
            "isotropic_limit_Q_Qd_zero": isotropic_limit_zero,
            "contains_Qdot": contains_qd,
            "contains_Q_squared": contains_q_squared,
        },
        "numeric_replay_table": rows,
        "gatekeeper_checks": {
            "p1980_eh_lapse_witness_present": p1980_pass,
            "r2_lapse_el_matches_symbolic_expected": el_matches_expected,
            "isotropic_limit_zero": isotropic_limit_zero,
            "r2_lapse_el_contains_Qdot": contains_qd,
            "r2_lapse_el_contains_Q_squared": contains_q_squared,
            "q_shear_positive_definite": q_positive_definite,
            "all_numeric_replay_samples_match_expected": all_samples_match,
            "r2_lapse_obligation_passed": theorem_pass,
        },
        "result_kind": "PASS_R2_LAPSE_VARIATION_OBLIGATION_EXPORTED" if theorem_pass else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "R2_LAPSE_EQUATION_HAS_NONTRIVIAL_QDOT_AND_Q2_SHEAR_TERMS",
            "RICCI2_RIEMANN2_GB_ADM_VARIATIONS_STILL_MISSING",
            "SPATIAL_EOM_TRANSPORT_STILL_MISSING",
            "BACKGROUND_INDEPENDENCE_REMAINS_OPEN",
            "GLOBAL_TOE_CLOSURE_STILL_OPEN",
        ],
        "theorem_export": {
            "positive_statement": "For trace-free diagonal Bianchi-I, R_BI-R_FRW=2*Q_shear/N^2.  The lapse Euler-Lagrange operator for the R^2 density difference is exactly -12*V*(6*H^2*Q - 2*H*Qdot + 4*Hdot*Q + Q^2)/N^4.",
            "obligation_statement": "Unlike the EH lapse term, the R^2 correction is not merely a constant multiple of -Q_shear; it contains Qdot and Q_shear^2.  Closure therefore requires the remaining curvature-squared sectors and spatial equations to be varied and combined with strict coefficients.",
            "not_licensed": [
                "Ricci^2 ADM variation",
                "Riemann^2 ADM variation",
                "Gauss-Bonnet ADM variation",
                "spatial Bianchi-I EOM cancellation",
                "global background-independence closure",
                "PO2/PO3 closure",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1981 is an R^2 lapse-equation obligation, not background-independence closure.  Its Qdot and Q^2 terms must not be silently dropped or absorbed without a full strict coefficient-level curvature-squared combination theorem.",
        "next_honest_step": "Compute the Ricci^2, Riemann^2, and Gauss-Bonnet ADM/Bianchi-I lapse Euler operators with N(t) retained, then combine all curvature-squared shear terms using the strict coefficients from P1972/P1853.",
        "lay_explanation": "Po dobrym wyniku dla zwykłej grawitacji Einsteina sprawdziliśmy pierwszy trudniejszy składnik: R^2.  On także reaguje na anizotropię, ale nie daje tylko prostego -Q_shear.  Pojawiają się dodatkowe terminy z pochodną ścinania i kwadratem ścinania.  To znaczy, że pełna teoria wymaga policzenia całego pakietu składników kwadratowych, a nie wolno zamknąć dowodu na samym Einsteinie-Hilbercie.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
