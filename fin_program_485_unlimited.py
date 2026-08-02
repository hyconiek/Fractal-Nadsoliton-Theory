#!/usr/bin/env python3
"""FIN P485: unbounded exact causal-axis ideal-membership computation.

There is intentionally no application-level time, CPU, or memory limit.
Operating-system and hardware limits still apply.  The script writes durable
stage markers before entering the Groebner computation.
"""

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
STATUS_PATH = ROOT / "FIN_Program_485_Unlimited_Status.json"
INPUT_PATH = ROOT / "FIN_Program_485_Exact_Ideal_Input.json"
RESULT_PATH = ROOT / "FIN_Program_485_Results.json"
BASIS_PATH = ROOT / "FIN_Program_485_Groebner_Basis.txt.gz"
REMAINDER_PATH = ROOT / "FIN_Program_485_Remainders.txt.gz"


def stamp(stage: str, **extra: object) -> None:
    payload = {
        "program": "P485",
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


def write_gzip(path: Path, lines: list[str]) -> str:
    data = ("\n\n".join(lines)+"\n").encode("utf-8")
    with gzip.open(path, "wb", compresslevel=6) as handle:
        handle.write(data)
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    started = time.monotonic()
    stamp("constructing_exact_reduced_branch")
    context = algebra.exact_context()
    reduced = algebra.reduced_branch_system(context)
    tangent = algebra.tangent_consistency(context)
    branch_substitution = reduced["branch_substitution"]
    targets = [
        sp.expand(value.subs(branch_substitution))
        for value in tangent["consistency_polynomials"]
    ]
    localizer = sp.symbols("pivot_inverse", real=True)
    variables = reduced["variables"][:-1] + (localizer, reduced["variables"][-1])
    generators = (
        list(reduced["equations"])
        + [reduced["alpha_polynomial"]]
        + [sp.expand(localizer*reduced["pivot_determinant"]-1)]
    )
    input_payload = {
        "program": "P485",
        "question": (
            "Do the exact positive standard-branch and localized active Riccati "
            "equations force all five causal-axis tangent consistency minors to zero?"
        ),
        "variables": [str(value) for value in variables],
        "active_pivot_rows": list(reduced["pivot"]["rows"]),
        "rho_reference_orbit": tangent["reference"],
        "generator_stats": algebra.polynomial_stats(generators, variables),
        "target_stats": algebra.polynomial_stats(targets, variables),
        "localization": "pivot_inverse*pivot_determinant-1",
        "branch": "125*(L-b)-36*alpha=0",
        "method": "SymPy exact QQ Groebner, grevlex, f5b",
        "limits": "none at application level",
    }
    INPUT_PATH.write_text(
        json.dumps(algebra.json_ready(input_payload), indent=2)+"\n",
        encoding="utf-8",
    )
    stamp(
        "groebner_running",
        variables=len(variables),
        generators=len(generators),
        targets=len(targets),
    )
    basis = sp.groebner(
        generators, *variables, order="grevlex", domain=sp.QQ, method="f5b"
    )
    basis_lines = [str(value.as_expr()) for value in basis.polys]
    basis_hash = write_gzip(BASIS_PATH, basis_lines)
    stamp("reducing_causal_axis_minors", basis_length=len(basis.polys))
    remainders = [sp.expand(basis.reduce(target)[1]) for target in targets]
    remainder_lines = [str(value) for value in remainders]
    remainder_hash = write_gzip(REMAINDER_PATH, remainder_lines)
    all_zero = all(value == 0 for value in remainders)
    p486 = json.loads((ROOT / "FIN_Program_486_Results.json").read_text(
        encoding="utf-8"
    ))
    premises_paid = bool(p486.get("orientation_premises_paid"))
    exact_phase_face = all_zero and premises_paid
    result = {
        "metadata": {
            "program": "P485",
            "execution_mode": "local unbounded exact algebra",
            "network_used": False,
            "laboratory_data_used": False,
            "external_audit_used": False,
        },
        "status": (
            "[Computer-assisted proof] all five causal-axis consistency minors "
            "belong to the localized exact positive-branch ideal; together with "
            "P486 this proves a nonzero exact causal tangent and refutes full-cone "
            "optimizer uniqueness"
            if exact_phase_face else
            "[Open after exact ideal test] the computed localized basis did not "
            "establish every required causal-axis consistency minor"
        ),
        "exact_phase_face_proved": exact_phase_face,
        "orientation_and_nonzero_reference_premises_paid_by_P486": premises_paid,
        "all_five_consistency_remainders_zero": all_zero,
        "remainder_zero_vector": [value == 0 for value in remainders],
        "remainder_stats": algebra.polynomial_stats(remainders, variables),
        "basis_length": len(basis.polys),
        "basis_is_zero_dimensional": basis.is_zero_dimensional,
        "basis_uncompressed_sha256": basis_hash,
        "remainders_uncompressed_sha256": remainder_hash,
        "elapsed_seconds": time.monotonic()-started,
        "theorem_if_passed": (
            "The exact P473 root lies on the positive standard branch and the "
            "localized pivot is nonzero. Vanishing of all consistency minors and "
            "the P486 nonzero reference coefficient give one exact rho and nonzero "
            "homogeneous-causal Q with X3 Q X3+K Q K/4=0. The P473 strict positivity "
            "margin then gives N+i*t*Q>0 for sufficiently small nonzero real t. "
            "Linearity preserves Riccati contact and objective value."
        ),
        "boundary": (
            "Even a successful result is a theorem about the declared finite comb. "
            "It is not a FIN selector, physical phase, unit, apparatus, or laboratory result."
        ),
        "new_object": "O189 Localized Exact Causal-Axis Ideal Certificate",
    }
    RESULT_PATH.write_text(
        json.dumps(algebra.json_ready(result), indent=2, ensure_ascii=False)+"\n",
        encoding="utf-8",
    )
    stamp("completed", result=str(RESULT_PATH), exact_phase_face=exact_phase_face)
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
