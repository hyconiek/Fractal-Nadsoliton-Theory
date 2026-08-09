#!/usr/bin/env python3
"""FIN P475: original full exact lexicographic elimination without limits."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import time
import traceback

import sympy as sp

import fin_phase_exact_algebra as algebra


ROOT = Path(__file__).resolve().parent
STATUS_PATH = ROOT / "FIN_Program_475_Unlimited_Status.json"
RESULT_PATH = ROOT / "FIN_Program_475_Unlimited_Results.json"
BASIS_PATH = ROOT / "FIN_Program_475_Unlimited_Lex_Basis.txt.gz"


def stamp(stage: str, **extra: object) -> None:
    payload = {
        "program": "P475-unlimited",
        "pid": os.getpid(),
        "stage": stage,
        "unix_time": time.time(),
        "no_application_time_limit": True,
        "no_application_memory_limit": True,
        **extra,
    }
    STATUS_PATH.write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")


def main() -> None:
    started = time.monotonic()
    stamp("constructing_full_14_variable_system")
    system = algebra.polynomial_system()
    alpha = sp.symbols("alpha", real=True)
    s1, s2, s3 = system["sines"]
    substitution = {
        s1: alpha,
        s2: 1-2*alpha**2,
        s3: 3*alpha-4*alpha**3,
    }
    A, B, u, L, a, b, c, d, e, f, g, h, i = system["variables"]
    variables = (A, B, u, a, b, c, d, e, f, g, h, i, alpha, L)
    generators = [sp.expand(value.subs(substitution)) for value in system["equations"]]
    generators.append(8*alpha**4-8*alpha**2+1)
    stamp("full_lex_groebner_running", variables=14, generators=14)
    basis = sp.groebner(
        generators, *variables, order="lex", domain=sp.QQ, method="f5b"
    )
    basis_text = "\n\n".join(str(value.as_expr()) for value in basis.polys)+"\n"
    raw = basis_text.encode("utf-8")
    with gzip.open(BASIS_PATH, "wb", compresslevel=6) as handle:
        handle.write(raw)
    univariate = [
        sp.Poly(value.as_expr(), L, domain=sp.QQ)
        for value in basis.polys
        if value.as_expr() != 0 and value.as_expr().free_symbols <= {L}
    ]
    result = {
        "metadata": {
            "program": "P475-unlimited",
            "execution_mode": "local unbounded full exact lexicographic elimination",
            "network_used": False,
            "laboratory_data_used": False,
            "external_audit_used": False,
        },
        "status": (
            "[Computer-assisted exact result] full elimination returned an L relation"
            if univariate else
            "[Completed without L relation] full exact basis has no L-only polynomial"
        ),
        "basis_length": len(basis.polys),
        "basis_is_zero_dimensional": basis.is_zero_dimensional,
        "basis_uncompressed_sha256": hashlib.sha256(raw).hexdigest(),
        "univariate_relations": [str(value.as_expr()) for value in univariate],
        "univariate_degrees": [int(value.degree()) for value in univariate],
        "elapsed_seconds": time.monotonic()-started,
        "boundary": (
            "A returned polynomial still requires branch isolation, irreducibility, "
            "and exact verification before it can be called the minimal polynomial of L."
        ),
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    stamp("completed", result=str(RESULT_PATH), relations=len(univariate))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    try:
        main()
    except BaseException as error:
        stamp(
            "failed",
            error_type=type(error).__name__,
            error=str(error),
            traceback=traceback.format_exc(),
        )
        raise
