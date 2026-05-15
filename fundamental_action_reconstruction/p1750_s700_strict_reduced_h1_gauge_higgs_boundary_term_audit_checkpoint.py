#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1749 = GEN / "p1749_s699_strict_reduced_h1_gauge_higgs_cross_variation_checkpoint.json"
OUT = GEN / "p1750_s700_strict_reduced_h1_gauge_higgs_boundary_term_audit_checkpoint.json"


def main() -> None:
    p1749 = json.loads(IN1749.read_text(encoding="utf-8"))

    x = sp.symbols("x", real=True)
    g = sp.Symbol("g", real=True)
    h = sp.Function("h")(x)

    raw_diff = sp.sympify(p1749["cross_variation"]["difference"])

    # Boundary-aware classification: check if obstruction is an exact derivative.
    # If yes, strict local equality fails, but weak-form equality can hold under boundary control.
    exact_derivative_candidate = sp.simplify(g * sp.diff(h, x))
    potential_primitive = sp.simplify(g * h)
    is_exact_derivative = sp.simplify(raw_diff - sp.diff(potential_primitive, x)) == 0

    boundary_clause = {
        "statement": "Integral of H1 difference vanishes if boundary term [g*h]_{∂Ω} = 0 is imposed.",
        "primitive": str(potential_primitive),
        "derivative": str(sp.diff(potential_primitive, x)),
        "exact_derivative_confirmed": bool(is_exact_derivative),
    }

    if is_exact_derivative:
        verdict = "BOUNDARY_SENSITIVE_OBSTRUCTION"
        strict_local_status = "FAIL_LOCAL_STRICT_ZERO"
        weak_form_status = "PASS_WEAK_FORM_WITH_BOUNDARY_CLAUSE"
    else:
        verdict = "OBSTRUCTION"
        strict_local_status = "FAIL_LOCAL_STRICT_ZERO"
        weak_form_status = "NOT_RECOVERED"

    payload = {
        "checkpoint": "P1750_S700",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "scope": "REDUCED_PROXY_ONLY_NOT_NONPROXY",
        "input_anchor": "p1749_s699",
        "h1_difference_raw": str(raw_diff),
        "boundary_term_audit": boundary_clause,
        "classification": {
            "verdict": verdict,
            "strict_local_status": strict_local_status,
            "weak_form_status": weak_form_status,
            "pass_zero_issued": False,
        },
        "policy_note": "No PASS_ZERO issued: weak-form salvage with boundary clause is not full strict-core closure.",
        "remaining_open": [
            "nonproxy_covariant_H1_A_mu_H_with_boundary_control",
            "explicit_covariant_E_A_mu_expression_nonproxy",
            "explicit_covariant_E_H_expression_nonproxy",
            "renormalization_closure_counterterm_stability",
            "cutkosky_unitarity_full_sector",
            "background_family_independence",
        ],
        "next_honest_step": "Uruchomić analogiczny boundary-aware H1 audit na jawnych nonproxy E_A^μ i E_H (4D) oraz dołączyć formalną kontrolę wyrazów brzegowych w kontrakcie P1732.",
        "lay_summary": "Wykryta przeszkoda nie jest losowa: to dokładna pochodna, więc znika po całkowaniu przy odpowiednim warunku brzegowym. To poprawia diagnozę, ale nadal nie zamyka pełnej teorii.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
