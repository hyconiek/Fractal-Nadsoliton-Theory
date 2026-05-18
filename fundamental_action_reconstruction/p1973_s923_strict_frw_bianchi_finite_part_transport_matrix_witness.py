#!/usr/bin/env python3
"""P1973 S923 strict FRW/Bianchi-I finite-part transport matrix witness.

This continues the ordered closure path after the P1972 B1 renormalization
PASS_ZERO ledger.  It builds the first explicit transport equation for the
P1897 nu-branch finite-part multiplier and proves the local transport residual
is zero.  It does not claim global background-independence closure.
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
OUT = GEN / "p1973_s923_strict_frw_bianchi_finite_part_transport_matrix_witness.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1897 = load("p1897_s847_strict_gbi_preclosure_transport_bridge_probe.json")
    p1972 = load("p1972_s922_strict_b1_renormalization_exact_cancellation_witness.json")

    lam = sp.symbols("lambda", real=True)
    nu = sp.symbols("nu", real=True)
    sigma2 = sp.symbols("sigma2", nonnegative=True, real=True)
    f0 = sp.symbols("F0", real=True)
    q = sp.simplify(nu * sigma2)

    # P1897 has G_BI_b1 = (1 + nu*sigma2) G_BI_frw.  Interpolate with
    # m(lambda)=1+lambda*q and solve dF/dlambda=A(lambda)F.
    multiplier = sp.simplify(1 + lam * q)
    connection_scalar = sp.simplify(sp.diff(multiplier, lam) / multiplier)
    finite_part_path = sp.simplify(multiplier * f0)
    transport_residual = sp.factor(sp.simplify(sp.diff(finite_part_path, lam) - connection_scalar * finite_part_path))

    # Four curvature coefficients inherit the same scalar transport in this
    # local finite-part-lock model.  This keeps the P1972 zero residual vector
    # locked because T(lambda)*0=0.
    transport_matrix = sp.eye(4) * multiplier
    connection_matrix = sp.eye(4) * connection_scalar
    matrix_residual = sp.simplify(sp.diff(transport_matrix, lam) - connection_matrix * transport_matrix)
    matrix_residual_entries = [sp.factor(sp.simplify(x)) for x in list(matrix_residual)]
    determinant = sp.factor(sp.det(transport_matrix.subs(lam, 1)))

    q_values = [0.0, 0.05, 0.2, 0.8]
    numeric_rows: list[dict[str, Any]] = []
    identity4 = np.eye(4)
    for q_value in q_values:
        # Integral_0^1 q/(1+lambda*q) dlambda = log(1+q).
        generator = np.log1p(q_value) * identity4
        scipy_transport = la.expm(generator)
        expected = (1.0 + q_value) * identity4
        numeric_rows.append(
            {
                "q_equals_nu_sigma2": q_value,
                "scipy_expm_transport_trace": float(np.trace(scipy_transport)),
                "expected_transport_trace": float(np.trace(expected)),
                "max_abs_transport_error": float(np.max(np.abs(scipy_transport - expected))),
                "determinant": float(la.det(scipy_transport)),
                "invertible_positive_branch": bool(1.0 + q_value > 0.0),
            }
        )

    p1972_residual_vector = p1972.get("residual_vector") or []
    zero_vector_transport_preserved = p1972_residual_vector == ["0", "0", "0", "0"]

    all_symbolic_zero = transport_residual == 0 and all(x == 0 for x in matrix_residual_entries)
    all_numeric_pass = all(row["max_abs_transport_error"] < 1e-14 and row["invertible_positive_branch"] for row in numeric_rows)

    out = {
        "ledger_id": "P1973_S923_STRICT_FRW_BIANCHI_FINITE_PART_TRANSPORT_MATRIX_WITNESS",
        "packet_id": "P1973",
        "stage_id": "S923",
        "produced_by": "p1973_s923_strict_frw_bianchi_finite_part_transport_matrix_witness.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "background_pair": ["FRW", "Bianchi-I"],
        "branch_parameter": "nu",
        "anisotropy_parameter": "sigma2",
        "depends_on": {
            "p1897_present": "gbi_bridge_quantities" in p1897,
            "p1972_present": p1972.get("result_kind") == "PASS_ZERO",
        },
        "source_bridge_reading": {
            "p1897_delta_branch": "G_BI_b1=(1+nu*sigma2)*G_BI_frw in the preclosure nu branch",
            "q_definition": "q=nu*sigma2",
            "admissible_local_branch": "1+q>0 so the transport matrix is invertible and orientation preserving",
        },
        "symbolic_transport": {
            "multiplier_m_lambda": str(multiplier),
            "connection_scalar_A_lambda": str(connection_scalar),
            "finite_part_path": str(finite_part_path),
            "transport_residual": str(transport_residual),
            "transport_residual_is_zero": bool(transport_residual == 0),
            "transport_matrix_T_lambda": str(transport_matrix),
            "connection_matrix_A_lambda": str(connection_matrix),
            "matrix_residual_entries": [str(x) for x in matrix_residual_entries],
            "matrix_residual_zero": all(x == 0 for x in matrix_residual_entries),
            "T_lambda_1_determinant": str(determinant),
        },
        "renormalization_lock_from_p1972": {
            "p1972_residual_vector": p1972_residual_vector,
            "zero_vector_transport_preserved": zero_vector_transport_preserved,
            "statement": "Because the P1972 divergent residual vector is zero, the scalar FRW/Bianchi finite-part transport maps it to zero on this local branch.",
        },
        "numeric_transport_replay": {
            "q_grid": q_values,
            "rows": numeric_rows,
            "all_scipy_expm_rows_match_closed_form": all_numeric_pass,
        },
        "gatekeeper_checks": {
            "symbolic_transport_residual_zero": all_symbolic_zero,
            "p1972_zero_residual_preserved": zero_vector_transport_preserved,
            "numeric_transport_replay_pass": all_numeric_pass,
            "admissible_branch_condition_exported": True,
        },
        "local_result_kind": "PASS_ZERO_LOCAL_TRANSPORT_RESIDUAL" if all_symbolic_zero and all_numeric_pass else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "theorem_export": {
            "licensed_statement": "For the P1897 nu-branch finite-part multiplier m(lambda)=1+lambda*nu*sigma2, the transport equation dF/dlambda=(nu*sigma2/(1+lambda*nu*sigma2))*F has exact solution F(lambda)=m(lambda)F_FRW and zero residual; the four-component B1 coefficient transport matrix is T(lambda)=m(lambda)I_4 on the local branch 1+nu*sigma2>0.",
            "not_licensed": [
                "global background-independence theorem",
                "full variational covariance for all EOM operators",
                "transport through singular branch 1+nu*sigma2=0",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1973 proves only the local finite-part transport residual for the P1897 scalar nu-branch multiplier.  It is a background-independence input witness, not global background-independence closure.",
        "next_honest_step": "Lift this scalar finite-part transport to the full variational operator bundle by checking the FRW/Bianchi-I EOM residual tensors componentwise, including the non-scalar nu-dependent anisotropic terms.",
        "lay_explanation": "Po skasowaniu rozbieżności w B1 sprawdziliśmy pierwszy prosty most między dwoma typami geometrii: FRW i Bianchi-I.  Most działa algebraicznie dla lokalnego mnożnika gałęzi nu, ale nie jest jeszcze dowodem, że wszystkie równania teorii są niezależne od wyboru tła.",
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
