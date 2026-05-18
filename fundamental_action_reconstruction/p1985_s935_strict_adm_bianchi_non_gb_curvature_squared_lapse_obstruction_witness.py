#!/usr/bin/env python3
"""P1985 S935 strict ADM/Bianchi-I non-GB curvature-squared lapse obstruction witness.

Combines the already exported ADM/Bianchi-I lapse operator differences for
R^2, Ricci^2 and Riemann^2 with strict B1 coefficients (P1853/P1972) and tests
whether the non-Gauss-Bonnet weighted combination cancels.

Scope: lapse operator structure only; strict lane only; no legacy-role transfer.
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
OUT = GEN / "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"


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
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")

    coeffs = (p1853.get("b1_symbolic_evaluation", {}).get("evaluated_coefficients", {}))
    a_r2 = sp.sympify(coeffs.get("a_R2", {}).get("symbolic", "0"))
    a_ric2 = sp.sympify(coeffs.get("a_Ric2", {}).get("symbolic", "0"))
    a_riem2 = sp.sympify(coeffs.get("a_Riem2", {}).get("symbolic", "0"))

    parse_locals = {
        name: sp.Symbol(name, real=True)
        for name in [
            "N", "Nd", "Ndd", "V", "H", "Hd", "Hdd",
            "sigma1", "sigma2", "dsigma1", "dsigma2", "d2sigma1", "d2sigma2",
            "Q", "Qd", "Qdd",
        ]
    }
    e_r2 = sp.sympify((p1981.get("r2_lapse_euler_operator", {}).get("EL_N_difference", "0")), locals=parse_locals)
    e_ric2 = sp.sympify((p1982.get("ricci2_lapse_euler_operator", {}).get("EL_N_difference", "0")), locals=parse_locals)
    e_riem2 = sp.sympify((p1983.get("riemann2_lapse_euler_operator", {}).get("EL_N_difference", "0")), locals=parse_locals)

    weighted_lapse_raw = sp.factor(sp.simplify(a_r2 * e_r2 + a_ric2 * e_ric2 + a_riem2 * e_riem2))

    s1, s2, sd1, sd2, sdd1, sdd2 = [parse_locals[k] for k in ("sigma1", "sigma2", "dsigma1", "dsigma2", "d2sigma1", "d2sigma2")]
    N, Nd, Ndd, V = [parse_locals[k] for k in ("N", "Nd", "Ndd", "V")]
    H, Hd, Hdd = [parse_locals[k] for k in ("H", "Hd", "Hdd")]
    q_sym, qd_sym = parse_locals["Q"], parse_locals["Qd"]
    q_expr = sp.factor(s1**2 + s1 * s2 + s2**2)
    qd_expr = sp.factor(2 * s1 * sd1 + s1 * sd2 + s2 * sd1 + 2 * s2 * sd2)

    weighted_lapse = sp.factor(sp.simplify(weighted_lapse_raw.subs({q_sym: q_expr, qd_sym: qd_expr})))

    anisotropic_nonzero = sp.simplify(weighted_lapse) != 0
    isotropic_zero = sp.simplify(weighted_lapse.subs({s1: 0, s2: 0, sd1: 0, sd2: 0, sdd1: 0, sdd2: 0})) == 0
    contains_high_derivatives = bool(weighted_lapse.has(sdd1) or weighted_lapse.has(sdd2) or weighted_lapse.has(Ndd))

    sample_points = [
        {N: sp.Rational(1, 1), Nd: sp.Rational(1, 20), Ndd: sp.Rational(-1, 200), V: sp.Rational(1, 1), H: sp.Rational(1, 1), Hd: sp.Rational(1, 10), Hdd: sp.Rational(1, 40), s1: sp.Rational(1, 10), s2: sp.Rational(-1, 20), sd1: sp.Rational(1, 100), sd2: sp.Rational(-1, 200), sdd1: sp.Rational(1, 1000), sdd2: sp.Rational(-1, 2000)},
        {N: sp.Rational(6, 5), Nd: sp.Rational(-1, 25), Ndd: sp.Rational(1, 200), V: sp.Rational(9, 8), H: sp.Rational(5, 4), Hd: sp.Rational(-1, 16), Hdd: sp.Rational(1, 64), s1: sp.Rational(3, 20), s2: sp.Rational(1, 20), sd1: sp.Rational(-1, 60), sd2: sp.Rational(1, 120), sdd1: sp.Rational(1, 840), sdd2: sp.Rational(-1, 1260)},
    ]
    replay = []
    vals: list[float] = []
    for i, point in enumerate(sample_points):
        val = sp.simplify(weighted_lapse.subs(point))
        replay.append({"sample_id": f"non_gb_lapse_sample_{i}", "weighted_lapse_EL": str(val), "is_zero": bool(val == 0)})
        vals.append(float(sp.N(val, 40)))

    nonzero_numeric_norm = float(la.norm(np.array(vals, dtype=float), ord=2))

    prereq = all([
        p1981.get("result_kind") == "PASS_R2_LAPSE_VARIATION_OBLIGATION_EXPORTED",
        p1982.get("result_kind") == "PASS_RICCI2_LAPSE_VARIATION_OBLIGATION_EXPORTED",
        p1983.get("result_kind") == "PASS_RIEMANN2_LAPSE_VARIATION_OBLIGATION_EXPORTED",
    ])
    witness_pass = prereq and anisotropic_nonzero and isotropic_zero and contains_high_derivatives and nonzero_numeric_norm > 0.0

    out = {
        "ledger_id": "P1985_S935_STRICT_ADM_BIANCHI_NON_GB_CURVATURE_SQUARED_LAPSE_OBSTRUCTION_WITNESS",
        "packet_id": "P1985",
        "stage_id": "S935",
        "produced_by": "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "strict_kernel": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
        "depends_on": {
            "p1981_present": p1981.get("result_kind") == "PASS_R2_LAPSE_VARIATION_OBLIGATION_EXPORTED",
            "p1982_present": p1982.get("result_kind") == "PASS_RICCI2_LAPSE_VARIATION_OBLIGATION_EXPORTED",
            "p1983_present": p1983.get("result_kind") == "PASS_RIEMANN2_LAPSE_VARIATION_OBLIGATION_EXPORTED",
            "p1853_coefficients_present": bool(coeffs),
        },
        "strict_coefficients": {"a_R2": str(a_r2), "a_Ric2": str(a_ric2), "a_Riem2": str(a_riem2)},
        "input_lapse_operators": {"EL_R2": str(e_r2), "EL_Ric2": str(e_ric2), "EL_Riem2": str(e_riem2)},
        "strict_symbolic_identifications": {
            "Q": str(q_expr),
            "Qd": str(qd_expr),
            "note": "P1981 symbols Q,Qd are identified with Bianchi-I shear invariants before combination."
        },
        "weighted_non_gb_lapse_operator": {
            "definition": "EL_nonGB := a_R2*EL_R2 + a_Ric2*EL_Ric2 + a_Riem2*EL_Riem2",
            "symbolic": str(weighted_lapse),
            "is_identically_zero": not anisotropic_nonzero,
            "isotropic_limit_zero": isotropic_zero,
            "contains_high_derivatives": contains_high_derivatives,
        },
        "numeric_replay_table": replay,
        "numeric_replay_l2_norm": nonzero_numeric_norm,
        "gatekeeper_checks": {
            "prerequisites_present": prereq,
            "non_gb_operator_nonzero": anisotropic_nonzero,
            "isotropic_limit_zero": isotropic_zero,
            "high_derivatives_present": contains_high_derivatives,
            "numeric_replay_nonzero_norm": nonzero_numeric_norm > 0.0,
            "obstruction_witness_passed": witness_pass,
        },
        "result_kind": "PASS_NON_GB_LAPSE_OBSTRUCTION_WITNESS" if witness_pass else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "STRICT_NON_GB_CURVATURE_SQUARED_LAPSE_RESIDUAL_NONZERO",
            "ANISOTROPIC_HIGH_DERIVATIVE_TERMS_PERSIST",
            "GLOBAL_COUNTERTERM_BRIDGE_STILL_OPEN",
        ],
        "false_pass_guard": "P1985 exports only the non-GB lapse residual witness; it does not close spatial EOM, background-independence, PO2/PO3, Cutkosky, or QW-2191.",
        "next_honest_step": "Project this strict non-GB lapse residual onto full Bianchi-I spatial EOM and test whether any admissible strict selector/provider class can cancel it without violating QW-2191 or kernel-split guardrails.",
        "lay_explanation": "Po zważeniu składników R^2, Ricci^2 i Riemann^2 prawdziwymi współczynnikami strict zostaje niezerowa reszta w równaniu lapse dla anizotropii. To znaczy, że sam trik Gaussa-Bonneta nie wystarcza i przeszkoda jest realna.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__("scipy").__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1985] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
