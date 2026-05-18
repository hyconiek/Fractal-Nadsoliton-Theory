#!/usr/bin/env python3
"""P1977 S927 strict positive-energy anisotropic provider bounded no-go.

This is the bounded theorem requested by the P1976 next step.  It proves a
conditional no-go under the current P1846/P1907 term-basis audit plus a
positive-energy requirement: no positive-energy anisotropic provider can be the
minimal componentwise canceller of the P1974 residual, because P1975 forces
rho_required = -Q_shear and Q_shear is positive definite for nonzero shear.

Scope: bounded to this route and these assumptions; not a global no-go theorem
for all future strict completions.
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
OUT = GEN / "p1977_s927_strict_positive_energy_anisotropic_provider_bounded_no_go.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1975 = load("p1975_s925_strict_minimal_anisotropic_source_obligation_and_energy_sign_audit.json")
    p1976 = load("p1976_s926_strict_ltotal_anisotropic_provider_nonexport_audit.json")

    s1, s2 = sp.symbols("sigma1 sigma2", real=True)
    q_shear = sp.factor(s1**2 + s1 * s2 + s2**2)
    rho_required = sp.sympify((p1975.get("source_obligation") or {}).get("rho_required", "0"), locals={"sigma1": s1, "sigma2": s2})
    rho_plus_q = sp.factor(sp.simplify(rho_required + q_shear))

    q_matrix = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    exact_eigs = sorted([sp.factor(ev) for ev in q_matrix.eigenvals().keys()], key=lambda x: float(x))
    scipy_eigs = la.eigvalsh(np.array(q_matrix.tolist(), dtype=float))
    q_positive_definite = all(float(ev) > 0 for ev in scipy_eigs)

    # Bounded contradiction: cancellation enforces rho=-Q; positive energy enforces rho>=0;
    # nonzero shear enforces Q>0, so rho<0 and rho>=0 cannot both hold.
    assumptions = {
        "A1_current_term_basis": "P1976 audits P1846/P1907 and exports no explicit anisotropic provider in current registries.",
        "A2_componentwise_cancellation": "P1975 minimal cancellation requires rho_provider = rho_required = -Q_shear.",
        "A3_positive_energy_provider": "Admissible matter/source provider must satisfy rho_provider >= 0 on physical branches.",
        "A4_nonzero_trace_free_shear": "(sigma1,sigma2) != (0,0), hence Q_shear > 0.",
    }
    contradiction_expression = "rho_provider = -Q_shear < 0 conflicts with rho_provider >= 0"

    sample_points = [
        {s1: sp.Rational(1, 10), s2: sp.Rational(-1, 20)},
        {s1: sp.Rational(1, 5), s2: sp.Rational(1, 10)},
        {s1: sp.Rational(-1, 8), s2: sp.Rational(1, 16)},
        {s1: sp.Rational(2, 7), s2: sp.Rational(-1, 9)},
    ]
    rows: list[dict[str, Any]] = []
    for idx, point in enumerate(sample_points):
        q_val = sp.simplify(q_shear.subs(point))
        rho_val = sp.simplify(rho_required.subs(point))
        rows.append(
            {
                "sample_id": f"bounded_no_go_sample_{idx}",
                "sigma1": str(point[s1]),
                "sigma2": str(point[s2]),
                "q_shear": str(q_val),
                "rho_required": str(rho_val),
                "positive_energy_compatible": bool(float(sp.N(rho_val, 30)) >= 0.0),
                "contradiction_detected": bool(float(sp.N(q_val, 30)) > 0.0 and float(sp.N(rho_val, 30)) < 0.0),
            }
        )

    p1976_nonexport = ((p1976.get("gatekeeper_checks") or {}).get("no_explicit_anisotropic_provider_in_current_registries") is True)
    all_samples_contradict = all(row["contradiction_detected"] for row in rows)
    theorem_pass = bool(rho_plus_q == 0 and q_positive_definite and p1976_nonexport and all_samples_contradict)

    out = {
        "ledger_id": "P1977_S927_STRICT_POSITIVE_ENERGY_ANISOTROPIC_PROVIDER_BOUNDED_NO_GO",
        "packet_id": "P1977",
        "stage_id": "S927",
        "produced_by": "p1977_s927_strict_positive_energy_anisotropic_provider_bounded_no_go.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1975_present": p1975.get("result_kind") == "OPEN_OBSTRUCTION_WITH_TRACE",
            "p1976_present": p1976.get("result_kind") == "OPEN_OBSTRUCTION_WITH_TRACE",
            "p1976_nonexport_audit_passed": p1976_nonexport,
        },
        "bounded_assumptions": assumptions,
        "symbolic_core": {
            "Q_shear": str(q_shear),
            "rho_required": str(rho_required),
            "rho_required_plus_Q_shear": str(rho_plus_q),
            "rho_required_equals_minus_Q_shear": bool(rho_plus_q == 0),
            "Q_matrix": str(q_matrix),
            "Q_eigenvalues_exact": [str(ev) for ev in exact_eigs],
            "Q_eigenvalues_scipy": [float(ev) for ev in scipy_eigs],
            "Q_positive_definite": q_positive_definite,
            "contradiction_expression": contradiction_expression,
        },
        "numeric_replay_table": rows,
        "gatekeeper_checks": {
            "p1976_current_basis_has_no_exported_provider": p1976_nonexport,
            "p1975_rho_equals_minus_q": bool(rho_plus_q == 0),
            "q_shear_positive_definite": q_positive_definite,
            "all_nonzero_shear_samples_violate_positive_energy": all_samples_contradict,
            "bounded_no_go_passed": theorem_pass,
        },
        "result_kind": "PASS_BOUNDED_NO_GO" if theorem_pass else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "POSITIVE_ENERGY_MINIMAL_ANISOTROPIC_PROVIDER_NO_GO_UNDER_CURRENT_BASIS",
            "BACKGROUND_INDEPENDENCE_REQUIRES_NONMINIMAL_TENSORIAL_OR_NEGATIVE_ENERGY_ROUTE",
            "GLOBAL_TOE_CLOSURE_STILL_OPEN",
        ],
        "theorem_export": {
            "bounded_no_go_statement": "Under the current P1846/P1907 term-basis audit, the P1975 minimal componentwise cancellation requirement, positive-energy rho_provider>=0, and nonzero trace-free shear, no positive-energy minimal anisotropic provider can cancel the P1974 residual, because rho_provider must equal -Q_shear<0.",
            "escape_routes": [
                "derive a non-minimal tensorial transport connection that avoids the minimal Delta_T obligation",
                "derive an explicitly gravitational/shear sector where the negative rho_required sign is allowed by the strict variational principle",
                "extend the strict L_total registry with a new theorem-grade anisotropic provider and rerun the sign audit",
            ],
            "not_licensed": [
                "global no-go theorem for all future strict completions",
                "background-independence closure",
                "negative-energy source admission",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1977 is a bounded no-go under explicit assumptions; it must not be read as a global impossibility theorem or as a closure of background-independence.",
        "next_honest_step": "Attempt the first escape route: derive a non-minimal tensorial transport connection for the Bianchi-I residual, or prove that the current strict variational operator bundle cannot provide one.",
        "lay_explanation": "Jeśli wymagamy zwykłej dodatniej energii, to brakujący anizotropowy składnik nie może być minimalnym źródłem naprawczym, bo matematyka wymaga energii ujemnej. To nie zabija teorii, ale wymusza trudniejszą drogę: nowe tensorowe przeniesienie albo specjalny sektor grawitacyjny z wyprowadzonym znakiem.",
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
