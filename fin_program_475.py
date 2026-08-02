#!/usr/bin/env python3
"""FIN Program P475: resource-bounded algebraic-value elimination audit."""

from __future__ import annotations

import csv
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

import numpy as np
import sympy as sp

import fin_programs_471_472_473 as p471


ROOT = Path(__file__).resolve().parent
RESULTS_PATH = ROOT / "FIN_Program_475_Results.json"
INVENTORY_PATH = ROOT / "FIN_Program_475_Algebraic_Inventory.csv"
WORKER_PATH = ROOT / "fin_program_475_elimination_worker.py"
WALL_SECONDS = 45


def json_ready(value: Any) -> Any:
    if isinstance(value, sp.Basic):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def exact_specialization() -> tuple[dict[str, Any], sp.Symbol, list[sp.Poly]]:
    system = p471.polynomial_system()
    alpha = sp.symbols("alpha", real=True)
    s1, s2, s3 = system["sines"]
    substitutions = {
        s1: alpha,
        s2: 1 - 2 * alpha**2,
        s3: 3 * alpha - 4 * alpha**3,
    }
    symbols = system["variables"] + (alpha,)
    polynomials = [
        sp.Poly(sp.expand(value.subs(substitutions)), *symbols, domain=sp.QQ)
        for value in system["equations"]
    ]
    return system, alpha, polynomials


def best_linear_pivot(
    system: dict[str, Any], alpha: sp.Symbol, polynomials: list[sp.Poly]
) -> dict[str, Any]:
    """Find the best-conditioned three-row local pivot for A,B,u.

    This is a scouting choice, not an exact elimination theorem.  Exact
    coefficient inspection separately proves every equation is affine in the
    three normalizer variables.
    """

    A, B, u = system["variables"][:3]
    root = np.load(ROOT / "FIN_Programs_471_472_473_P473_Root_Box.npz")["center"]
    substitutions = {
        symbol: float(value) for symbol, value in zip(system["variables"], root)
    }
    substitutions[alpha] = float(np.sin(np.pi / 8))
    expressions = [poly.as_expr() for poly in polynomials]
    coefficient_matrix = sp.Matrix([
        [sp.diff(value, variable) for variable in (A, B, u)]
        for value in expressions
    ])
    numeric = np.asarray(coefficient_matrix.subs(substitutions), dtype=float)
    candidates: list[tuple[float, tuple[int, int, int]]] = []
    for rows in itertools.combinations(range(13), 3):
        determinant = float(np.linalg.det(numeric[list(rows), :]))
        candidates.append((abs(determinant), rows))
    score, rows = max(candidates)
    exact_determinant = sp.Poly(
        sp.expand(coefficient_matrix[list(rows), :].det()),
        *system["variables"][3:], alpha,
        domain=sp.QQ,
    )
    return {
        "rows": rows,
        "floating_absolute_determinant": score,
        "exact_pivot_determinant_is_nonzero_polynomial": not exact_determinant.is_zero,
        "exact_pivot_determinant_total_degree": exact_determinant.total_degree(),
        "exact_pivot_determinant_term_count": len(exact_determinant.terms()),
    }


def bounded_worker() -> dict[str, Any]:
    started = time.monotonic()
    environment = dict(os.environ)
    environment.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    try:
        completed = subprocess.run(
            [sys.executable, str(WORKER_PATH)],
            cwd=ROOT,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=WALL_SECONDS,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        return {
            "completed": False,
            "termination": "wall_clock_timeout",
            "wall_limit_seconds": WALL_SECONDS,
            "elapsed_seconds": time.monotonic() - started,
            "stdout_tail": (error.stdout or "")[-1000:] if isinstance(error.stdout, str) else "",
            "stderr_tail": (error.stderr or "")[-1000:] if isinstance(error.stderr, str) else "",
        }
    elapsed = time.monotonic() - started
    parsed = None
    if completed.returncode == 0:
        try:
            parsed = json.loads(completed.stdout.strip().splitlines()[-1])
        except (json.JSONDecodeError, IndexError):
            parsed = None
    return {
        "completed": bool(parsed and parsed.get("completed")),
        "termination": "completed" if parsed else "worker_nonzero_or_unparseable",
        "return_code": completed.returncode,
        "wall_limit_seconds": WALL_SECONDS,
        "worker_cpu_limit_seconds": 40,
        "worker_address_space_limit_bytes": 3 * 1024**3,
        "elapsed_seconds": elapsed,
        "worker_result": parsed,
        "stdout_tail": completed.stdout[-1000:],
        "stderr_tail": completed.stderr[-1000:],
    }


def main() -> None:
    system, alpha, polynomials = exact_specialization()
    variables = system["variables"]
    normalizer_variables = set(variables[:3])
    rows = []
    affine_all = True
    for index, polynomial in enumerate(polynomials):
        expression = polynomial.as_expr()
        normalizer_degrees = {
            str(variable): sp.degree(expression, variable) for variable in variables[:3]
        }
        affine = all((degree or 0) <= 1 for degree in normalizer_degrees.values())
        no_products = all(
            sp.diff(expression, left, right) == 0
            for left in variables[:3]
            for right in variables[:3]
        )
        affine_all = affine_all and affine and no_products
        rows.append({
            "equation": index,
            "total_degree_Q_alpha": polynomial.total_degree(),
            "monomial_count": len(polynomial.terms()),
            "degree_alpha": sp.degree(expression, alpha),
            "affine_in_A_B_u": affine and no_products,
        })
    with INVENTORY_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    pivot = best_linear_pivot(system, alpha, polynomials)
    worker = bounded_worker()
    completed = worker["completed"]
    univariate = (
        worker.get("worker_result", {}).get("univariate_in_L", []) if completed else []
    )
    minimal_polynomial_obtained = bool(univariate)
    result = {
        "metadata": {
            "program": "P475",
            "execution_mode": "local analytical and resource-bounded computational research",
            "network_used": False,
            "laboratory_data_used": False,
            "external_audit_used": False,
            "hardware_envelope": "Intel i3-10110U, 16 GB RAM",
        },
        "question": "Does the exact O181 optimum L have a tractable algebraic minimal polynomial under the declared local envelope?",
        "status": (
            "[Computer-assisted result] lexicographic elimination completed and produced a univariate L relation"
            if minimal_polynomial_obtained else
            "[Open; resource-bounded no-go] no univariate minimal-polynomial candidate was obtained within the declared exact-elimination envelope"
        ),
        "exact_field_reduction": {
            "primitive_element": "alpha=sin(pi/8)",
            "primitive_polynomial": "8*alpha**4 - 8*alpha**2 + 1",
            "sine_reconstruction": {
                "sin(pi/8)": "alpha",
                "sin(pi/4)": "1 - 2*alpha**2",
                "sin(3*pi/8)": "3*alpha - 4*alpha**3",
            },
            "field_degree": 4,
            "exact": True,
        },
        "system_inventory": {
            "equations_before_field_relation": len(polynomials),
            "unknown_operator_coordinates": len(variables),
            "rational_elimination_variables_with_alpha": len(variables) + 1,
            "maximum_total_degree_after_Q_specialization": max(poly.total_degree() for poly in polynomials),
            "maximum_monomial_count": max(len(poly.terms()) for poly in polynomials),
            "total_monomial_occurrences": sum(len(poly.terms()) for poly in polynomials),
            "all_equations_affine_in_normalizer_coordinates_A_B_u": affine_all,
            "naive_cubic_Bezout_bound_over_Q_alpha": 3**13,
            "warning": "The Bezout number is only a crude projective complexity bound and is not a root count for the admitted positive branch.",
        },
        "local_linear_pivot": pivot,
        "bounded_lexicographic_attempt": worker,
        "minimal_polynomial_candidate_obtained": minimal_polynomial_obtained,
        "candidate_relations": univariate,
        "proved": [
            "The three trigonometric coefficients lie in the degree-four algebraic field Q(alpha).",
            "The exact Riccati system is polynomial over Q(alpha) and affine in A, B, and u.",
            "The declared bounded local SymPy lexicographic workflow either completed as recorded or was forcibly stopped by its stated resource envelope.",
        ],
        "not_proved": [
            "No claim is made that L lacks a minimal polynomial; L is expected to be algebraic if the certified isolated branch is a zero-dimensional algebraic component.",
            "No irreducibility, algebraic degree, global ideal dimension, or complete solution count is proved unless a univariate relation was actually returned.",
            "A timeout or resource kill is not a mathematical impossibility theorem and does not exclude modular, resultant, homotopy, or structure-aware elimination.",
        ],
        "falsification": {
            "alternative_checked": "Exact one-generator trigonometric specialization and an exact affine normalizer-variable audit were performed before lexicographic elimination.",
            "failure_mode": "A future structure-aware or modular reconstruction can refute the tractability boundary by producing and verifying a univariate relation inside a comparable envelope.",
        },
        "new_object": "O187 Resource-Bounded Algebraic Elimination Burden Certificate",
        "boundary": (
            "P475 concerns a dimensionless exact operator optimum only. It does not close the legacy-to-strict bridge, QW-2191, physical units, apparatus, or laboratory evidence."
        ),
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(result), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(json_ready(result), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
