#!/usr/bin/env python3
"""P1987 S937 strict non-GB residual term-classification witness.

Next honest step after P1986: classify the surviving non-GB lapse residual into
term families (H^2*Q, H*Qd, Hd*Q, Q^2, d2sigma*sigma, dsigma^2, Ndd*sigma^2,
Nd^2*sigma^2, ...), export exact coefficients, and verify the highest-risk
families are genuinely nonzero in strict lane.
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
OUT = GEN / "p1987_s937_strict_non_gb_residual_term_classification_witness.json"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1986 = load("p1986_s936_strict_adm_bianchi_non_gb_residual_decomposition_witness.json")
    p1985 = load("p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json")

    parse = {
        n: sp.Symbol(n, real=True)
        for n in [
            "N", "Nd", "Ndd", "V", "H", "Hd", "Hdd",
            "sigma1", "sigma2", "dsigma1", "dsigma2", "d2sigma1", "d2sigma2",
        ]
    }
    s1, s2, sd1, sd2, sdd1, sdd2 = [parse[k] for k in ("sigma1", "sigma2", "dsigma1", "dsigma2", "d2sigma1", "d2sigma2")]
    N, Nd, Ndd, V, H, Hd = [parse[k] for k in ("N", "Nd", "Ndd", "V", "H", "Hd")]

    expr = sp.sympify(p1986.get("non_gb_residual", {}).get("symbolic", "0"), locals=parse)
    if expr == 0:
        expr = sp.sympify(p1985.get("weighted_non_gb_lapse_operator", {}).get("symbolic", "0"), locals={**parse, "Q": sp.Symbol("Q"), "Qd": sp.Symbol("Qd")})
        Q = s1**2 + s1 * s2 + s2**2
        Qd = 2 * s1 * sd1 + s1 * sd2 + s2 * sd1 + 2 * s2 * sd2
        expr = sp.simplify(expr.subs({sp.Symbol("Q"): Q, sp.Symbol("Qd"): Qd}))

    scaled = sp.factor(sp.together(expr * N**6 / V))

    families = {
        "H2_Q": [H**2 * s1**2, H**2 * s1 * s2, H**2 * s2**2],
        "H_Qd": [H * s1 * sd1, H * s1 * sd2, H * s2 * sd1, H * s2 * sd2],
        "Hd_Q": [Hd * s1**2, Hd * s1 * s2, Hd * s2**2],
        "Q2": [s1**4, s1**3 * s2, s1**2 * s2**2, s1 * s2**3, s2**4],
        "d2sigma_sigma": [sdd1 * s1, sdd1 * s2, sdd2 * s1, sdd2 * s2],
        "dsigma_sq": [sd1**2, sd1 * sd2, sd2**2],
        "Ndd_sigma_sq": [Ndd * s1**2, Ndd * s1 * s2, Ndd * s2**2],
        "Nd2_sigma_sq": [Nd**2 * s1**2, Nd**2 * s1 * s2, Nd**2 * s2**2],
    }

    coeffs: dict[str, str] = {}
    nonzero_flags: dict[str, bool] = {}
    expanded = sp.expand(scaled)
    for key, monomials in families.items():
        parts = []
        nz = False
        for m in monomials:
            c = sp.factor(expanded.coeff(m))
            parts.append(f"coeff[{sp.sstr(m)}]={sp.sstr(c)}")
            nz = nz or bool(sp.simplify(c) != 0)
        coeffs[key] = "; ".join(parts)
        nonzero_flags[key] = nz

    # Robust fallback extraction directly by differentiating with dummy parameterization.
    x = sp.Symbol("x")
    check_expr = sp.expand(scaled.subs({
        s1: x * s1, s2: x * s2, sd1: x * sd1, sd2: x * sd2, sdd1: x * sdd1, sdd2: x * sdd2
    }))
    order2_present = bool(sp.expand(check_expr).coeff(x, 2) != 0)
    order3_present = bool(sp.expand(check_expr).coeff(x, 3) != 0)
    order4_present = bool(sp.expand(check_expr).coeff(x, 4) != 0)

    sample = {
        N: sp.Rational(1, 1), Nd: sp.Rational(1, 20), Ndd: sp.Rational(-1, 200),
        V: sp.Rational(1, 1), H: sp.Rational(1, 1), Hd: sp.Rational(1, 10),
        s1: sp.Rational(1, 10), s2: sp.Rational(-1, 20), sd1: sp.Rational(1, 100), sd2: sp.Rational(-1, 200),
        sdd1: sp.Rational(1, 1000), sdd2: sp.Rational(-1, 2000),
    }
    val = sp.simplify(expr.subs(sample))
    l2 = float(la.norm(np.array([float(sp.N(val, 40))], dtype=float), ord=2))

    key_families = ["H2_Q", "H_Qd", "Hd_Q", "Q2", "d2sigma_sigma", "Ndd_sigma_sq"]
    key_nonzero = all(nonzero_flags[k] for k in key_families)

    out = {
        "ledger_id": "P1987_S937_STRICT_NON_GB_RESIDUAL_TERM_CLASSIFICATION_WITNESS",
        "packet_id": "P1987",
        "stage_id": "S937",
        "produced_by": Path(__file__).name,
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1986_present": p1986.get("result_kind") == "PASS_NON_GB_DECOMPOSITION_OBSTRUCTION_WITNESS",
            "p1985_present": p1985.get("result_kind") == "PASS_NON_GB_LAPSE_OBSTRUCTION_WITNESS",
        },
        "residual_expression": str(expr),
        "scaled_expression_N6_over_V": str(scaled),
        "term_family_coefficients": coeffs,
        "term_family_nonzero_flags": nonzero_flags,
        "anisotropy_order_presence": {
            "order2_present": order2_present,
            "order3_present": order3_present,
            "order4_present": order4_present,
        },
        "numeric_probe": {
            "sample": {str(k): str(v) for k, v in sample.items()},
            "residual_value": str(val),
            "l2_norm": l2,
        },
        "gatekeeper_checks": {
            "key_term_families_nonzero": key_nonzero,
            "anisotropy_orders_2_and_4_present": bool(order2_present and order4_present),
            "order3_present": order3_present,
            "numeric_probe_nonzero": l2 > 0.0,
        },
        "result_kind": "PASS_NON_GB_TERM_CLASSIFICATION_WITNESS" if (key_nonzero and order2_present and order4_present and l2 > 0.0) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "next_honest_step": "Project the classified term families onto full spatial Bianchi-I EOM equations and test if any strict admissible provider can cancel each family without violating QW-2191.",
        "lay_explanation": "Rozpisaliśmy pozostały błąd na klasy składników. To pokazuje dokładnie, które typy wyrazów (np. z pochodnymi ścinania i lapse) blokują domknięcie i gdzie trzeba celować następnym dowodem.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1987] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
