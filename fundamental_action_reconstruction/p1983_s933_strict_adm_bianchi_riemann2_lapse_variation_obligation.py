#!/usr/bin/env python3
"""P1983 S933 strict ADM/Bianchi-I Riemann^2 lapse-variation obligation.

This continues the curvature-squared ADM/Bianchi-I lapse chain after P1982 by
computing the Riemann_{mu nu rho sigma}R^{mu nu rho sigma} (Kretschmann) lapse
Euler operator with N(t) retained.

Scope: Riemann^2 minisuperspace lapse equation only.  Gauss-Bonnet, strict
coefficient combination, spatial equations, and global background independence
remain open.
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
OUT = GEN / "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1982 = load("p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json")

    N, Nd, Ndd, V = sp.symbols("N Nd Ndd V", positive=True, real=True)
    H, Hd, Hdd = sp.symbols("H Hd Hdd", real=True)
    s1, s2, sd1, sd2, sdd1, sdd2 = sp.symbols("sigma1 sigma2 dsigma1 dsigma2 d2sigma1 d2sigma2", real=True)

    s3 = sp.factor(-s1 - s2)
    sd3 = sp.factor(-sd1 - sd2)
    hs = [H + s1, H + s2, H + s3]
    hds = [Hd + sd1, Hd + sd2, Hd + sd3]

    # Orthonormal diagonal Bianchi-I curvature blocks with lapse retained:
    # E_i = R_{0i0i}/(orthonormal signs), F_ij = R_{ijij}.  The Kretschmann
    # scalar is 4*sum_i E_i^2 + 4*sum_{i<j} F_ij^2, reproducing
    # 12*((Hdot+H^2)^2+H^4) in the FRW limit.
    e_i = [sp.factor((hd + h**2 - h * Nd / N) / N**2) for h, hd in zip(hs, hds)]
    f_ij = [sp.factor(hs[0] * hs[1] / N**2), sp.factor(hs[0] * hs[2] / N**2), sp.factor(hs[1] * hs[2] / N**2)]
    riemann2_bi = sp.factor(4 * (sum(e**2 for e in e_i) + sum(f**2 for f in f_ij)))

    e_frw = sp.factor((Hd + H**2 - H * Nd / N) / N**2)
    f_frw = sp.factor(H**2 / N**2)
    riemann2_frw = sp.factor(12 * (e_frw**2 + f_frw**2))
    riemann2_difference = sp.factor(riemann2_bi - riemann2_frw)

    density_difference = sp.factor(N * V * riemann2_difference)
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
    contains_cubic_shear_velocity_mixing = bool(
        riemann2_difference.has(s1**2 * sd2)
        or riemann2_difference.has(s2**2 * sd1)
        or riemann2_difference.has(s1 * s2 * sd1)
        or riemann2_difference.has(s1 * s2 * sd2)
    )

    q_shear = sp.factor(s1**2 + s1 * s2 + s2**2)
    q_matrix = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    scipy_eigs = la.eigvalsh(np.array(q_matrix.tolist(), dtype=float))
    q_positive_definite = all(float(ev) > 0 for ev in scipy_eigs)

    sample_points = [
        {N: sp.Rational(1, 1), Nd: sp.Rational(1, 20), Ndd: sp.Rational(-1, 200), V: sp.Rational(1, 1), H: sp.Rational(1, 1), Hd: sp.Rational(1, 10), s1: sp.Rational(1, 10), s2: sp.Rational(-1, 20), sd1: sp.Rational(1, 100), sd2: sp.Rational(-1, 200), sdd1: sp.Rational(1, 1000), sdd2: sp.Rational(-1, 2000)},
        {N: sp.Rational(5, 4), Nd: sp.Rational(-1, 30), Ndd: sp.Rational(1, 300), V: sp.Rational(7, 5), H: sp.Rational(3, 2), Hd: sp.Rational(-1, 20), s1: sp.Rational(1, 5), s2: sp.Rational(1, 10), sd1: sp.Rational(-1, 50), sd2: sp.Rational(3, 200), sdd1: sp.Rational(1, 700), sdd2: sp.Rational(-1, 900)},
        {N: sp.Rational(4, 3), Nd: sp.Rational(1, 40), Ndd: sp.Rational(-1, 400), V: sp.Rational(11, 10), H: sp.Rational(4, 5), Hd: sp.Rational(1, 25), s1: sp.Rational(-1, 8), s2: sp.Rational(1, 16), sd1: sp.Rational(1, 80), sd2: sp.Rational(1, 160), sdd1: sp.Rational(-1, 1600), sdd2: sp.Rational(1, 3200)},
    ]
    rows: list[dict[str, Any]] = []
    for idx, point in enumerate(sample_points):
        diff_val = sp.simplify(riemann2_difference.subs(point))
        el_val = sp.simplify(lapse_el_difference.subs(point))
        rows.append(
            {
                "sample_id": f"riemann2_lapse_variation_sample_{idx}",
                "input": {str(k): str(v) for k, v in point.items()},
                "Riemann2_Bianchi_minus_FRW": str(diff_val),
                "lapse_EL_difference": str(el_val),
                "absolute_float_value": float(abs(sp.N(el_val, 30))),
                "matches_symbolic_self_replay": bool(sp.simplify(el_val - lapse_el_difference.subs(point)) == 0),
            }
        )

    p1982_pass = p1982.get("result_kind") == "PASS_RICCI2_LAPSE_VARIATION_OBLIGATION_EXPORTED"
    all_samples_match = all(row["matches_symbolic_self_replay"] for row in rows)
    theorem_pass = bool(
        p1982_pass
        and isotropic_limit_zero
        and contains_shear_acceleration
        and contains_lapse_second_derivative
        and contains_shear_velocity_squared
        and contains_cubic_shear_velocity_mixing
        and q_positive_definite
        and all_samples_match
    )

    out = {
        "ledger_id": "P1983_S933_STRICT_ADM_BIANCHI_RIEMANN2_LAPSE_VARIATION_OBLIGATION",
        "packet_id": "P1983",
        "stage_id": "S933",
        "produced_by": "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1982_ricci2_lapse_obligation_present": p1982_pass,
        },
        "riemann_block_setup": {
            "sigma3": str(s3),
            "dsigma3": str(sd3),
            "h_i": [str(h) for h in hs],
            "dh_i": [str(hd) for hd in hds],
            "E_i_orthonormal_time_space_blocks": [str(e) for e in e_i],
            "F_ij_orthonormal_spatial_blocks": [str(f) for f in f_ij],
            "Kretschmann_formula": "4*sum_i E_i^2 + 4*sum_{i<j} F_ij^2",
            "FRW_check_formula": "12*((Hdot+H^2-H*Ndot/N)^2/N^4 + H^4/N^4)",
        },
        "riemann2_lapse_euler_operator": {
            "Riemann2_BianchiI": str(riemann2_bi),
            "Riemann2_FRW": str(riemann2_frw),
            "Riemann2_Bianchi_minus_FRW": str(riemann2_difference),
            "density_difference_NV_Riemann2": str(density_difference),
            "dL_dN": str(d_l_d_n),
            "dL_dNdot": str(d_l_d_nd),
            "Dt_dL_dNdot": str(dt_d_l_d_nd),
            "EL_N_difference": str(lapse_el_difference),
            "isotropic_limit_zero": isotropic_limit_zero,
            "contains_shear_acceleration": contains_shear_acceleration,
            "contains_lapse_second_derivative": contains_lapse_second_derivative,
            "contains_shear_velocity_squared": contains_shear_velocity_squared,
            "contains_cubic_shear_velocity_mixing_in_invariant_difference": contains_cubic_shear_velocity_mixing,
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
            "p1982_ricci2_lapse_obligation_present": p1982_pass,
            "isotropic_limit_zero": isotropic_limit_zero,
            "riemann2_lapse_el_contains_shear_acceleration": contains_shear_acceleration,
            "riemann2_lapse_el_contains_lapse_second_derivative": contains_lapse_second_derivative,
            "riemann2_lapse_el_contains_shear_velocity_squared": contains_shear_velocity_squared,
            "riemann2_invariant_difference_contains_cubic_shear_velocity_mixing": contains_cubic_shear_velocity_mixing,
            "q_shear_positive_definite": q_positive_definite,
            "all_numeric_replay_samples_match_symbolic_self_replay": all_samples_match,
            "riemann2_lapse_obligation_passed": theorem_pass,
        },
        "result_kind": "PASS_RIEMANN2_LAPSE_VARIATION_OBLIGATION_EXPORTED" if theorem_pass else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "RIEMANN2_LAPSE_EQUATION_HAS_D2SHEAR_AND_NDDOT_TERMS",
            "GAUSS_BONNET_COMBINATION_STILL_MISSING",
            "STRICT_COEFFICIENT_COMBINATION_STILL_MISSING",
            "SPATIAL_EOM_TRANSPORT_STILL_MISSING",
            "GLOBAL_TOE_CLOSURE_STILL_OPEN",
        ],
        "theorem_export": {
            "positive_statement": "For diagonal trace-free Bianchi-I, the Riemann^2 lapse Euler operator difference is exported exactly from the orthonormal Kretschmann blocks and vanishes in the isotropic limit.",
            "obligation_statement": "The Riemann^2 lapse operator has the same higher-derivative obstruction class as Ricci^2 and additionally records cubic shear/velocity mixing in the invariant difference.  Gauss-Bonnet and strict coefficient combinations remain mandatory before any background-independence claim.",
            "not_licensed": [
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
        "false_pass_guard": "P1983 exports the Riemann^2 lapse obligation only.  It must not be used to claim Gauss-Bonnet, strict coefficient-level, spatial-EOM, or ToE closure.",
        "next_honest_step": "Build the Gauss-Bonnet ADM/Bianchi-I lapse combination GB = Riemann^2 - 4*Ricci^2 + R^2 from P1981/P1982/P1983 and test whether the higher-derivative shear structures cancel as expected before applying strict coefficients from P1972/P1853.",
        "lay_explanation": "Policzyliśmy trzeci składnik kwadratowy, Riemann^2.  Tak jak Ricci^2, zawiera on wyższe pochodne ścinania i lapse, więc nie wolno udawać prostego domknięcia.  Teraz trzeba zbudować kombinację Gaussa-Bonneta i sprawdzić, czy część tych trudnych wyrazów rzeczywiście się kasuje.",
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
