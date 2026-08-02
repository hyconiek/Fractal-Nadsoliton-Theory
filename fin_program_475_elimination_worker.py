#!/usr/bin/env python3
"""Bounded exact lexicographic elimination worker for FIN Program P475.

The parent process imposes a wall-clock timeout.  This worker additionally
sets conservative CPU and address-space ceilings before asking SymPy for a
Groebner basis.  A killed or timed-out worker is evidence about this declared
implementation/resource envelope only, never an algebraic impossibility.
"""

from __future__ import annotations

import json
import resource
import time

import sympy as sp

import fin_programs_471_472_473 as p471


CPU_SECONDS = 40
ADDRESS_SPACE_BYTES = 3 * 1024**3


def main() -> None:
    resource.setrlimit(resource.RLIMIT_CPU, (CPU_SECONDS, CPU_SECONDS + 1))
    resource.setrlimit(
        resource.RLIMIT_AS, (ADDRESS_SPACE_BYTES, ADDRESS_SPACE_BYTES)
    )
    system = p471.polynomial_system()
    alpha = sp.symbols("alpha", real=True)
    s1, s2, s3 = system["sines"]
    substitutions = {
        s1: alpha,
        s2: 1 - 2 * alpha**2,
        s3: 3 * alpha - 4 * alpha**3,
    }
    A, B, u, L, a, b, c, d, e, f, g, h, i = system["variables"]
    # Lexicographic order eliminates all entries and alpha before L.
    order = (A, B, u, a, b, c, d, e, f, g, h, i, alpha, L)
    polynomials = [sp.expand(value.subs(substitutions)) for value in system["equations"]]
    polynomials.append(8 * alpha**4 - 8 * alpha**2 + 1)
    started = time.monotonic()
    basis = sp.groebner(polynomials, *order, order="lex", domain=sp.QQ)
    elapsed = time.monotonic() - started
    univariate = [
        sp.Poly(value, L, domain=sp.QQ)
        for value in basis.polys
        if value.as_expr().free_symbols <= {L} and value.as_expr() != 0
    ]
    print(json.dumps({
        "completed": True,
        "elapsed_seconds": elapsed,
        "basis_length": len(basis.polys),
        "zero_dimensional": basis.is_zero_dimensional,
        "univariate_in_L": [str(value.as_expr()) for value in univariate],
        "univariate_degrees": [value.degree() for value in univariate],
        "maximum_resident_set_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }))


if __name__ == "__main__":
    main()
