#!/usr/bin/env python3
"""FIN P487: unbounded structure-aware exact elimination for the value L."""

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
STATUS_PATH = ROOT / "FIN_Program_487_Unlimited_Status.json"
INPUT_PATH = ROOT / "FIN_Program_487_Reduced_Ideal_Input.json"
RESULT_PATH = ROOT / "FIN_Program_487_Results.json"
BASIS_PATH = ROOT / "FIN_Program_487_Lex_Groebner_Basis.txt.gz"


def stamp(stage: str, **extra: object) -> None:
    payload = {
        "program": "P487",
        "pid": os.getpid(),
        "stage": stage,
        "unix_time": time.time(),
        "no_application_time_limit": True,
        "no_application_memory_limit": True,
        **extra,
    }
    STATUS_PATH.write_text(
        json.dumps(algebra.json_ready(payload), indent=2)+"\n", encoding="utf-8"
    )


def main() -> None:
    started = time.monotonic()
    stamp("constructing_structure_aware_reduction")
    context = algebra.exact_context()
    reduced = algebra.reduced_branch_system(context)
    localizer = sp.symbols("pivot_inverse", real=True)
    base_variables = reduced["variables"]
    variables = base_variables[:-1] + (localizer, base_variables[-1])
    generators = (
        list(reduced["equations"])
        + [reduced["alpha_polynomial"]]
        + [sp.expand(localizer*reduced["pivot_determinant"]-1)]
    )
    input_payload = {
        "program": "P487",
        "variables": [str(value) for value in variables],
        "active_pivot_rows": list(reduced["pivot"]["rows"]),
        "remaining_active_rows": reduced["remaining_rows"],
        "generator_stats": algebra.polynomial_stats(generators, variables),
        "normalizer_numerator_stats": algebra.polynomial_stats(
            reduced["normalizer_numerators"], variables
        ),
        "branch": "b=L-36*alpha/125",
        "localization": "pivot_inverse*pivot_determinant-1",
        "method": "SymPy exact QQ Groebner, lex, f5b; L ordered last",
        "limits": "none at application level",
    }
    INPUT_PATH.write_text(
        json.dumps(algebra.json_ready(input_payload), indent=2)+"\n",
        encoding="utf-8",
    )
    stamp("lex_groebner_running", variables=len(variables), generators=len(generators))
    basis = sp.groebner(
        generators, *variables, order="lex", domain=sp.QQ, method="f5b"
    )
    basis_text = ("\n\n".join(str(value.as_expr()) for value in basis.polys)+"\n")
    raw = basis_text.encode("utf-8")
    with gzip.open(BASIS_PATH, "wb", compresslevel=6) as handle:
        handle.write(raw)
    L = base_variables[-1]
    univariate = [
        sp.Poly(value.as_expr(), L, domain=sp.QQ)
        for value in basis.polys
        if value.as_expr() != 0 and value.as_expr().free_symbols <= {L}
    ]
    result = {
        "metadata": {
            "program": "P487",
            "execution_mode": "local unbounded structure-aware exact elimination",
            "network_used": False,
            "laboratory_data_used": False,
            "external_audit_used": False,
        },
        "status": (
            "[Computer-assisted exact result] the structure-aware localized "
            "positive-branch elimination returned at least one univariate relation in L"
            if univariate else
            "[Completed without L relation] the exact reduced lexicographic basis "
            "contains no nonzero polynomial solely in L"
        ),
        "basis_length": len(basis.polys),
        "basis_is_zero_dimensional": basis.is_zero_dimensional,
        "basis_uncompressed_sha256": hashlib.sha256(raw).hexdigest(),
        "univariate_relations": [str(value.as_expr()) for value in univariate],
        "univariate_degrees": [int(value.degree()) for value in univariate],
        "univariate_factorizations": [str(sp.factor(value.as_expr())) for value in univariate],
        "elapsed_seconds": time.monotonic()-started,
        "boundary": (
            "Any returned relation belongs to the localized algebraic branch. "
            "A minimal-polynomial claim additionally requires isolating the P473 "
            "factor, irreducibility, and exact substitution/interval verification."
        ),
        "new_object": "O190 Structure-Aware Localized O181 Elimination System",
    }
    RESULT_PATH.write_text(
        json.dumps(algebra.json_ready(result), indent=2, ensure_ascii=False)+"\n",
        encoding="utf-8",
    )
    stamp("completed", result=str(RESULT_PATH), relations=len(univariate))
    print(json.dumps(algebra.json_ready(result), indent=2, ensure_ascii=False))


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
