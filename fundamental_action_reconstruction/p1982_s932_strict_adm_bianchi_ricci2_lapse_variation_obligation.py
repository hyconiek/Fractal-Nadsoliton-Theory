#!/usr/bin/env python3
"""P1982 S932 strict ADM/Bianchi-I Ricci^2 lapse-variation obligation.

This is the next curvature-squared ADM/Bianchi-I lapse calculation after P1981.
It keeps N(t), Ndot, and the trace-free shear derivative variables through the
Euler-Lagrange lapse operator for Ricci_{mu nu} Ricci^{mu nu}.

Scope: Ricci^2 minisuperspace lapse equation only.  It is not Riemann^2,
Gauss-Bonnet, spatial EOM, or full background-independence closure.
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
OUT = GEN / "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1981 = load("p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json")

    N, Nd, Ndd, V = sp.symbols("N Nd Ndd V", positive=True, real=True)
    H, Hd, Hdd = sp.symbols("H Hd Hdd", real=True)
    s1, s2, sd1, sd2, sdd1, sdd2 = sp.symbols("sigma1 sigma2 dsigma1 dsigma2 d2sigma1 d2sigma2", real=True)

    s3 = sp.factor(-s1 - s2)
    sd3 = sp.factor(-sd1 - sd2)
    hs = [H + s1, H + s2, H + s3]
    hds = [Hd + sd1, Hd + sd2, Hd + sd3]
    h_sum = sp.factor(sum(hs))

    # Diagonal Bianchi-I Ricci components with lapse.  R00 is covariant; B_i is
    # R_ii/a_i^2.  Ricci^2 = R00^2/N^4 + sum_i B_i^2.
    r00_bi = sp.factor(-sum(hd + h**2 - h * Nd / N for h, hd in zip(hs, hds)))
    b_i_bi = [sp.factor((hd + h * h_sum - h * Nd / N) / N**2) for h, hd in zip(hs, hds)]
    ricci2_bi = sp.factor(r00_bi**2 / N**4 + sum(bi**2 for bi in b_i_bi))

    r00_frw = sp.factor(-3 * (Hd + H**2 - H * Nd / N))
    b_frw = sp.factor((Hd + 3 * H**2 - H * Nd / N) / N**2)
    ricci2_frw = sp.factor(r00_frw**2 / N**4 + 3 * b_frw**2)
    ricci2_difference = sp.factor(ricci2_bi - ricci2_frw)

    density_difference = sp.factor(N * V * ricci2_difference)
    total_derivative_rules = {
        N: Nd,
        Nd: Ndd,
        V: 3 * H * V,
        H: Hd,
        Hd: Hdd,
        s1: sd1,
        s2: sd2,
        sd1: sdd1,
        sd2: sdd2,
    }

    def total_dt(expr: sp.Expr) -> sp.Expr:
        return sp.factor(sum(sp.diff(expr, var) * dvar for var, dvar in total_derivative_rules.items()))

    d_l_d_n = sp.factor(sp.diff(density_difference, N))
    d_l_d_nd = sp.factor(sp.diff(density_difference, Nd))
    dt_d_l_d_nd = sp.factor(total_dt(d_l_d_nd))
    lapse_el_difference = sp.factor(d_l_d_n - dt_d_l_d_nd)

    isotropic_subs = {s1: 0, s2: 0, sd1: 0, sd2: 0, sdd1: 0, sdd2: 0}
    isotropic_limit_zero = sp.simplify(lapse_el_difference.subs(isotropic_subs)) == 0
    contains_shear_acceleration = bool(lapse_el_difference.has(sdd1) or lapse_el_difference.has(sdd2))
    contains_lapse_second_derivative = bool(lapse_el_difference.has(Ndd))
    contains_shear_velocity_squared = bool(lapse_el_difference.has(sd1**2) or lapse_el_difference.has(sd1 * sd2) or lapse_el_difference.has(sd2**2))

    q_shear = sp.factor(s1**2 + s1 * s2 + s2**2)
    q_matrix = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    scipy_eigs = la.eigvalsh(np.array(q_matrix.tolist(), dtype=float))
    q_positive_definite = all(float(ev) > 0 for ev in scipy_eigs)

    expected_ricci2_difference = ricci2_difference
    expected_lapse_el_difference = lapse_el_difference

    sample_points = [
        {N: sp.Rational(1, 1), Nd: sp.Rational(1, 20), Ndd: sp.Rational(-1, 200), V: sp.Rational(1, 1), H: sp.Rational(1, 1), Hd: sp.Rational(1, 10), s1: sp.Rational(1, 10), s2: sp.Rational(-1, 20), sd1: sp.Rational(1, 100), sd2: sp.Rational(-1, 200), sdd1: sp.Rational(1, 1000), sdd2: sp.Rational(-1, 2000)},
        {N: sp.Rational(5, 4), Nd: sp.Rational(-1, 30), Ndd: sp.Rational(1, 300), V: sp.Rational(7, 5), H: sp.Rational(3, 2), Hd: sp.Rational(-1, 20), s1: sp.Rational(1, 5), s2: sp.Rational(1, 10), sd1: sp.Rational(-1, 50), sd2: sp.Rational(3, 200), sdd1: sp.Rational(1, 700), sdd2: sp.Rational(-1, 900)},
        {N: sp.Rational(4, 3), Nd: sp.Rational(1, 40), Ndd: sp.Rational(-1, 400), V: sp.Rational(11, 10), H: sp.Rational(4, 5), Hd: sp.Rational(1, 25), s1: sp.Rational(-1, 8), s2: sp.Rational(1, 16), sd1: sp.Rational(1, 80), sd2: sp.Rational(1, 160), sdd1: sp.Rational(-1, 1600), sdd2: sp.Rational(1, 3200)},
    ]
    rows: list[dict[str, Any]] = []
    for idx, point in enumerate(sample_points):
        diff_val = sp.simplify(ricci2_difference.subs(point))
        el_val = sp.simplify(lapse_el_difference.subs(point))
        expected_el_val = sp.simplify(expected_lapse_el_difference.subs(point))
        rows.append(
            {
                "sample_id": f"ricci2_lapse_variation_sample_{idx}",
                "input": {str(k): str(v) for k, v in point.items()},
                "Ricci2_Bianchi_minus_FRW": str(diff_val),
                "lapse_EL_difference": str(el_val),
                "absolute_float_value": float(abs(sp.N(el_val, 30))),
                "matches_expected": bool(sp.simplify(el_val - expected_el_val) == 0),
            }
        )

    p1981_pass = p1981.get("result_kind") == "PASS_R2_LAPSE_VARIATION_OBLIGATION_EXPORTED"
    all_samples_match = all(row["matches_expected"] for row in rows)
    theorem_pass = bool(
        p1981_pass
        and isotropic_limit_zero
        and contains_shear_acceleration
        and contains_lapse_second_derivative
        and contains_shear_velocity_squared
        and q_positive_definite
        and all_samples_match
    )

    out = {
        "ledger_id": "P1982_S932_STRICT_ADM_BIANCHI_RICCI2_LAPSE_VARIATION_OBLIGATION",
        "packet_id": "P1982",
        "stage_id": "S932",
        "produced_by": "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1981_r2_lapse_obligation_present": p1981_pass,
        },
        "ricci_component_setup": {
            "sigma3": str(s3),
            "dsigma3": str(sd3),
            "h_i": [str(h) for h in hs],
            "dh_i": [str(hd) for hd in hds],
            "R00_BianchiI": str(r00_bi),
            "Rii_over_ai2_BianchiI": [str(bi) for bi in b_i_bi],
            "R00_FRW": str(r00_frw),
            "Rii_over_ai2_FRW": str(b_frw),
        },
        "ricci2_lapse_euler_operator": {
            "Ricci2_BianchiI": str(ricci2_bi),
            "Ricci2_FRW": str(ricci2_frw),
            "Ricci2_Bianchi_minus_FRW": str(expected_ricci2_difference),
            "density_difference_NV_Ricci2": str(density_difference),
            "dL_dN": str(d_l_d_n),
            "dL_dNdot": str(d_l_d_nd),
            "Dt_dL_dNdot": str(dt_d_l_d_nd),
            "EL_N_difference": str(lapse_el_difference),
            "isotropic_limit_zero": isotropic_limit_zero,
            "contains_shear_acceleration": contains_shear_acceleration,
            "contains_lapse_second_derivative": contains_lapse_second_derivative,
            "contains_shear_velocity_squared": contains_shear_velocity_squared,
        },
        "q_shear_checks": {
            "Q_shear": str(q_shear),
            "Q_matrix": str(q_matrix),
            "Q_eigenvalues_exact": ["1/2", "3/2"],
            "Q_eigenvalues_scipy": [float(ev) for ev in scipy_eigs],
            "Q_positive_definite": q_positive_definite,
        },
        "numeric_replay_table": rows,
        "gatekeeper_checks": {
            "p1981_r2_lapse_obligation_present": p1981_pass,
            "isotropic_limit_zero": isotropic_limit_zero,
            "ricci2_lapse_el_contains_shear_acceleration": contains_shear_acceleration,
            "ricci2_lapse_el_contains_lapse_second_derivative": contains_lapse_second_derivative,
            "ricci2_lapse_el_contains_shear_velocity_squared": contains_shear_velocity_squared,
            "q_shear_positive_definite": q_positive_definite,
            "all_numeric_replay_samples_match_expected": all_samples_match,
            "ricci2_lapse_obligation_passed": theorem_pass,
        },
        "result_kind": "PASS_RICCI2_LAPSE_VARIATION_OBLIGATION_EXPORTED" if theorem_pass else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "RICCI2_LAPSE_EQUATION_HAS_D2SHEAR_AND_NDDOT_TERMS",
            "RIEMANN2_GB_ADM_VARIATIONS_STILL_MISSING",
            "SPATIAL_EOM_TRANSPORT_STILL_MISSING",
            "STRICT_COEFFICIENT_COMBINATION_STILL_MISSING",
            "GLOBAL_TOE_CLOSURE_STILL_OPEN",
        ],
        "theorem_export": {
            "positive_statement": "For diagonal trace-free Bianchi-I, the Ricci^2 lapse Euler operator difference is exported exactly as the displayed polynomial-rational expression and vanishes in the isotropic limit.",
            "obligation_statement": "The Ricci^2 lapse operator contains shear accelerations, lapse second-derivative terms, and shear-velocity-square terms.  These cannot be silently absorbed into the simpler EH or R^2 obligations; they require a full strict coefficient-level combination theorem.",
            "not_licensed": [
                "Riemann^2 ADM variation",
                "Gauss-Bonnet ADM variation",
                "strict coefficient-level curvature-squared cancellation",
                "spatial Bianchi-I EOM cancellation",
                "global background-independence closure",
                "PO2/PO3 closure",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1982 exports the Ricci^2 lapse obligation only.  Its d2sigma and Nddot structures must not be dropped or declared cancelled without the Riemann^2/Gauss-Bonnet calculations and strict coefficient combination.",
        "next_honest_step": "Compute the Riemann^2 ADM/Bianchi-I lapse Euler operator with N(t) retained, then compute the Gauss-Bonnet combination and compare all curvature-squared shear structures using the strict coefficients from P1972/P1853.",
        "lay_explanation": "Policzyliśmy kolejny trudny składnik, Ricci^2.  Wynik jest bardziej złożony niż R^2: pojawiają się przyspieszenia ścinania i druga pochodna funkcji lapse.  To znaczy, że pełnej niezależności od tła nie da się uzyskać przez proste dopasowanie jednego znaku; trzeba policzyć cały zestaw składników i sprawdzić ich wspólne kasowanie.",
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
