#!/usr/bin/env python3
"""P1989 S939 strict minimal provider family no-go witness.

Next honest step after P1988: test a minimal strict anisotropic-provider ansatz
family and verify it cannot cancel all outside-EH non-GB channels simultaneously.

Important scope: this is only a no-go for a minimal provider class; not a global
impossibility theorem for all provider classes.
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
OUT = GEN / "p1989_s939_strict_minimal_provider_family_no_go_witness.json"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1987 = load("p1987_s937_strict_non_gb_residual_term_classification_witness.json")
    p1988 = load("p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json")

    s1, s2, sd1, sd2, sdd1, sdd2 = sp.symbols("sigma1 sigma2 dsigma1 dsigma2 d2sigma1 d2sigma2", real=True)
    N, Nd, Ndd, V, H, Hd = sp.symbols("N Nd Ndd V H Hd", real=True)

    locals_map = {
        "N": N, "Nd": Nd, "Ndd": Ndd, "V": V, "H": H, "Hd": Hd,
        "sigma1": s1, "sigma2": s2, "dsigma1": sd1, "dsigma2": sd2, "d2sigma1": sdd1, "d2sigma2": sdd2,
        "pi": sp.pi, "log": sp.log, "ln": sp.log,
    }
    residual = sp.sympify(p1987.get("scaled_expression_N6_over_V", "0"), locals=locals_map)

    # Minimal strict provider class: EH-like anisotropic channels only
    # P = c1*(H*sigma linear) + c2*(dsigma linear) + c3*Q.
    c1, c2, c3 = sp.symbols("c1 c2 c3", real=True)
    lin_sigma = 3 * H * s1 + 3 * H * s2
    lin_dsigma = sd1 + sd2
    Q = s1**2 + s1 * s2 + s2**2
    provider = sp.expand(c1 * lin_sigma + c2 * lin_dsigma + c3 * Q)

    target = sp.expand(residual - provider)

    # Coefficients of outside-EH channels that must vanish for closure.
    channel_coeffs = {
        "d2sigma1_sigma1": sp.expand(target).coeff(sdd1 * s1),
        "d2sigma2_sigma2": sp.expand(target).coeff(sdd2 * s2),
        "Ndd_sigma1_sq": sp.expand(target).coeff(Ndd * s1**2),
        "Q_quartic_s1_4": sp.expand(target).coeff(s1**4),
    }

    # Solve exact system channel_coeffs == 0 for (c1,c2,c3).
    eqs = [sp.Eq(v, 0) for v in channel_coeffs.values()]
    sol = sp.solve(eqs, (c1, c2, c3), dict=True)

    no_solution_minimal = len(sol) == 0

    # Numeric replay with random rational substitutions for symbols and c's set to best-fit 0.
    sample = {
        N: sp.Rational(1, 1), Nd: sp.Rational(1, 20), Ndd: sp.Rational(-1, 200),
        V: sp.Rational(1, 1), H: sp.Rational(1, 1), Hd: sp.Rational(1, 10),
        s1: sp.Rational(1, 10), s2: sp.Rational(-1, 20), sd1: sp.Rational(1, 100), sd2: sp.Rational(-1, 200),
        sdd1: sp.Rational(1, 1000), sdd2: sp.Rational(-1, 2000), c1: 0, c2: 0, c3: 0,
    }
    val = sp.simplify(target.subs(sample))
    l2 = float(la.norm(np.array([float(sp.N(val, 40))], dtype=float), ord=2))

    gate = {
        "p1987_present": p1987.get("result_kind") == "PASS_NON_GB_TERM_CLASSIFICATION_WITNESS",
        "p1988_present": p1988.get("result_kind") == "PASS_NON_GB_SPATIAL_PROJECTION_OBSTRUCTION_WITNESS",
        "outside_eh_channels_reported": len((p1988.get("family_projection") or {}).get("outside_eh_family_capacity", [])) > 0,
        "minimal_provider_has_no_exact_solution": no_solution_minimal,
        "numeric_residual_nonzero": l2 > 0.0,
    }

    out = {
        "ledger_id": "P1989_S939_STRICT_MINIMAL_PROVIDER_FAMILY_NO_GO_WITNESS",
        "packet_id": "P1989",
        "stage_id": "S939",
        "produced_by": Path(__file__).name,
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1987_present": gate["p1987_present"],
            "p1988_present": gate["p1988_present"],
        },
        "minimal_provider_family": {
            "definition": "P_min = c1*(3H(sigma1+sigma2)) + c2*(dsigma1+dsigma2) + c3*Q",
            "Q": str(Q),
            "unknowns": ["c1", "c2", "c3"],
        },
        "outside_channel_equations": {k: str(v) for k, v in channel_coeffs.items()},
        "symbolic_solver": {
            "equations": [str(e) for e in eqs],
            "solutions_for_c1_c2_c3": sol,
            "no_solution_in_minimal_family": no_solution_minimal,
        },
        "numeric_probe": {
            "sample": {str(k): str(v) for k, v in sample.items()},
            "residual_value": str(val),
            "l2_norm": l2,
        },
        "gatekeeper_checks": gate,
        "result_kind": "PASS_MINIMAL_PROVIDER_NO_GO_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "MINIMAL_EH_LIKE_PROVIDER_INSUFFICIENT",
            "OUTSIDE_EH_CHANNELS_PERSIST",
            "STRICT_SELECTOR_AUGMENTED_PROVIDER_CLASS_STILL_OPEN",
        ],
        "false_pass_guard": "P1989 is only a no-go witness for one minimal EH-like provider family; it does not exclude broader strict admissible provider classes.",
        "next_honest_step": "Introduce one extended strict provider ansatz class with explicit selector-premise label (non-strict unless bridged) and re-run channel-equation solvability family-by-family.",
        "lay_explanation": "Sprawdziliśmy najprostszy legalny kandydat źródła anizotropowego. Okazało się, że nie potrafi on skasować wszystkich trudnych składników residualu. To znaczy, że potrzeba bogatszego (ale nadal rygorystycznego) mechanizmu.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1989] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
