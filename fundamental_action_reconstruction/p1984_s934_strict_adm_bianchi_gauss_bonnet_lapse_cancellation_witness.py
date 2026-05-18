#!/usr/bin/env python3
"""P1984 S934 strict ADM/Bianchi-I Gauss-Bonnet lapse cancellation witness.

Builds the Gauss-Bonnet combination from the already executed R^2, Ricci^2 and
Riemann^2 ADM/Bianchi-I lapse chain and verifies the lapse Euler operator
cancels exactly in the minisuperspace witness.

Scope: lapse equation for the GB combination only.  This is not a full
coefficient-level curvature-squared cancellation theorem, not a spatial EOM
closure, and not ToE closure.
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
OUT = GEN / "p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1981 = load("p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json")
    p1982 = load("p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json")
    p1983 = load("p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json")

    N, Nd, Ndd, V = sp.symbols("N Nd Ndd V", positive=True, real=True)
    H, Hd, Hdd = sp.symbols("H Hd Hdd", real=True)
    s1, s2, sd1, sd2, sdd1, sdd2 = sp.symbols("sigma1 sigma2 dsigma1 dsigma2 d2sigma1 d2sigma2", real=True)

    s3 = sp.factor(-s1 - s2)
    sd3 = sp.factor(-sd1 - sd2)
    hs = [H + s1, H + s2, H + s3]
    hds = [Hd + sd1, Hd + sd2, Hd + sd3]
    h_sum = sp.factor(sum(hs))
    q_shear = sp.factor(s1**2 + s1 * s2 + s2**2)

    # R^2 difference from the lapse-retaining Ricci scalar.
    a_frw = sp.factor(3 * Hd + 6 * H**2 - 3 * H * Nd / N)
    r_frw = sp.factor(2 * a_frw / N**2)
    r_bi = sp.factor(2 * (a_frw + q_shear) / N**2)
    r2_diff = sp.factor(r_bi**2 - r_frw**2)

    # Ricci^2 difference.
    r00_bi = sp.factor(-sum(hd + h**2 - h * Nd / N for h, hd in zip(hs, hds)))
    b_i_bi = [sp.factor((hd + h * h_sum - h * Nd / N) / N**2) for h, hd in zip(hs, hds)]
    ricci2_bi = sp.factor(r00_bi**2 / N**4 + sum(bi**2 for bi in b_i_bi))
    r00_frw = sp.factor(-3 * (Hd + H**2 - H * Nd / N))
    b_frw = sp.factor((Hd + 3 * H**2 - H * Nd / N) / N**2)
    ricci2_frw = sp.factor(r00_frw**2 / N**4 + 3 * b_frw**2)
    ricci2_diff = sp.factor(ricci2_bi - ricci2_frw)

    # Riemann^2/Kretschmann difference.
    e_i = [sp.factor((hd + h**2 - h * Nd / N) / N**2) for h, hd in zip(hs, hds)]
    f_ij = [sp.factor(hs[0] * hs[1] / N**2), sp.factor(hs[0] * hs[2] / N**2), sp.factor(hs[1] * hs[2] / N**2)]
    riemann2_bi = sp.factor(4 * (sum(e**2 for e in e_i) + sum(f**2 for f in f_ij)))
    e_frw = sp.factor((Hd + H**2 - H * Nd / N) / N**2)
    f_frw = sp.factor(H**2 / N**2)
    riemann2_frw = sp.factor(12 * (e_frw**2 + f_frw**2))
    riemann2_diff = sp.factor(riemann2_bi - riemann2_frw)

    gb_diff = sp.factor(riemann2_diff - 4 * ricci2_diff + r2_diff)
    gb_density_diff = sp.factor(N * V * gb_diff)

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

    d_l_d_n = sp.factor(sp.diff(gb_density_diff, N))
    d_l_d_nd = sp.factor(sp.diff(gb_density_diff, Nd))
    dt_d_l_d_nd = sp.factor(total_dt(d_l_d_nd))
    gb_lapse_el = sp.factor(d_l_d_n - dt_d_l_d_nd)

    higher_derivative_cancelled = not (gb_lapse_el.has(sdd1) or gb_lapse_el.has(sdd2) or gb_lapse_el.has(Ndd))
    lapse_el_zero = sp.simplify(gb_lapse_el) == 0
    isotropic_limit_zero = sp.simplify(gb_lapse_el.subs({s1: 0, s2: 0, sd1: 0, sd2: 0, sdd1: 0, sdd2: 0})) == 0
    gb_density_nonzero_polynomial = sp.simplify(gb_diff) != 0

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
        gb_density_val = sp.simplify(gb_density_diff.subs(point))
        gb_el_val = sp.simplify(gb_lapse_el.subs(point))
        rows.append(
            {
                "sample_id": f"gauss_bonnet_lapse_sample_{idx}",
                "input": {str(k): str(v) for k, v in point.items()},
                "GB_density_difference": str(gb_density_val),
                "GB_lapse_EL": str(gb_el_val),
                "GB_lapse_EL_is_zero": bool(gb_el_val == 0),
            }
        )

    p1981_pass = p1981.get("result_kind") == "PASS_R2_LAPSE_VARIATION_OBLIGATION_EXPORTED"
    p1982_pass = p1982.get("result_kind") == "PASS_RICCI2_LAPSE_VARIATION_OBLIGATION_EXPORTED"
    p1983_pass = p1983.get("result_kind") == "PASS_RIEMANN2_LAPSE_VARIATION_OBLIGATION_EXPORTED"
    all_samples_zero = all(row["GB_lapse_EL_is_zero"] for row in rows)
    theorem_pass = bool(
        p1981_pass
        and p1982_pass
        and p1983_pass
        and lapse_el_zero
        and isotropic_limit_zero
        and higher_derivative_cancelled
        and gb_density_nonzero_polynomial
        and q_positive_definite
        and all_samples_zero
    )

    out = {
        "ledger_id": "P1984_S934_STRICT_ADM_BIANCHI_GAUSS_BONNET_LAPSE_CANCELLATION_WITNESS",
        "packet_id": "P1984",
        "stage_id": "S934",
        "produced_by": "p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1981_r2_lapse_obligation_present": p1981_pass,
            "p1982_ricci2_lapse_obligation_present": p1982_pass,
            "p1983_riemann2_lapse_obligation_present": p1983_pass,
        },
        "gauss_bonnet_combination": {
            "definition": "G_GB = Riemann^2 - 4*Ricci^2 + R^2",
            "R2_difference": str(r2_diff),
            "Ricci2_difference": str(ricci2_diff),
            "Riemann2_difference": str(riemann2_diff),
            "GB_difference": str(gb_diff),
            "GB_density_difference_NV": str(gb_density_diff),
            "GB_density_nonzero_polynomial_before_variation": gb_density_nonzero_polynomial,
        },
        "gb_lapse_euler_operator": {
            "dL_dN": str(d_l_d_n),
            "dL_dNdot": str(d_l_d_nd),
            "Dt_dL_dNdot": str(dt_d_l_d_nd),
            "EL_N_GB_difference": str(gb_lapse_el),
            "EL_N_GB_difference_is_zero": lapse_el_zero,
            "isotropic_limit_zero": isotropic_limit_zero,
            "higher_derivative_terms_cancelled_in_lapse_EL": higher_derivative_cancelled,
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
            "all_input_lapse_obligations_present": bool(p1981_pass and p1982_pass and p1983_pass),
            "gb_lapse_euler_operator_zero": lapse_el_zero,
            "higher_derivative_shear_lapse_terms_cancel_in_gb_lapse_EL": higher_derivative_cancelled,
            "gb_density_difference_nonzero_before_variation": gb_density_nonzero_polynomial,
            "all_numeric_replay_samples_zero": all_samples_zero,
            "q_shear_positive_definite": q_positive_definite,
            "gb_lapse_cancellation_witness_passed": theorem_pass,
        },
        "result_kind": "PASS_GB_LAPSE_CANCELLATION_WITNESS" if theorem_pass else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "GAUSS_BONNET_LAPSE_EL_CANCELS_IN_MINISUPERSPACE",
            "STRICT_NON_GB_COEFFICIENT_COMBINATION_STILL_MISSING",
            "SPATIAL_EOM_TRANSPORT_STILL_MISSING",
            "GLOBAL_BACKGROUND_INDEPENDENCE_STILL_OPEN",
            "GLOBAL_TOE_CLOSURE_STILL_OPEN",
        ],
        "theorem_export": {
            "positive_statement": "In the diagonal trace-free Bianchi-I ADM minisuperspace witness, the Gauss-Bonnet combination Riemann^2 - 4*Ricci^2 + R^2 has identically zero lapse Euler operator, although the density difference is not the zero polynomial before variation.",
            "interpretation": "The higher-derivative shear/lapse structures exposed separately by P1981/P1982/P1983 cancel in the GB lapse Euler equation, matching the expected topological behavior at the minisuperspace lapse level.",
            "not_licensed": [
                "strict coefficient-level cancellation for non-GB R^2/Ricci^2/Riemann^2 terms",
                "spatial Bianchi-I EOM cancellation",
                "global background-independence closure",
                "PO2/PO3 closure",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1984 proves only the Gauss-Bonnet lapse Euler cancellation in this ADM/Bianchi-I minisuperspace witness.  It must not be upgraded to full curvature-squared, spatial-EOM, or ToE closure.",
        "next_honest_step": "Use P1981/P1982/P1983/P1984 to form the strict non-GB curvature-squared lapse operator with coefficients a_R2, a_Ric2, and a_Riem2 from P1972/P1853, then test whether the remaining shear structures cancel or persist.",
        "lay_explanation": "Zbudowaliśmy kombinację Gaussa-Bonneta z trzech policzonych klocków.  Trudne wyrazy z wyższymi pochodnymi kasują się w równaniu lapse, co jest dobrym znakiem i zgodne z topologicznym charakterem Gaussa-Bonneta.  To jednak nie zamyka teorii, bo pozostałe nie-Gauss-Bonnetowe kombinacje ze ścisłymi współczynnikami nadal trzeba sprawdzić.",
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
