#!/usr/bin/env python3
"""Build the serialized exact certificate consumed by the P480 checker.

This builder may use the scientific Python stack.  The emitted checker and
certificate are designed so that replay itself uses only the Python standard
library.
"""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

import fin_programs_471_472_473 as p471


ROOT = Path(__file__).resolve().parent
CERTIFICATE = ROOT / "FIN_Program_480_Standalone_Certificate.json"


def serialize_fraction(value: Fraction | sp.Rational) -> str:
    return str(Fraction(int(value.p), int(value.q))) if isinstance(value, sp.Rational) else str(value)


def serialize_interval(value: tuple[Fraction, Fraction]) -> list[str]:
    return [str(value[0]), str(value[1])]


def serialize_polynomial(polynomial: sp.Poly) -> list[dict[str, Any]]:
    return [
        {
            "exponents": list(exponents),
            "coefficient": serialize_fraction(coefficient),
        }
        for exponents, coefficient in polynomial.terms()
    ]


def serialize_matrix(matrix: Any) -> list[list[str]]:
    """Serialize either a nested sequence or a SymPy Matrix by rows."""

    rows = matrix.tolist() if isinstance(matrix, sp.MatrixBase) else matrix
    return [[str(value) for value in row] for row in rows]


def main() -> None:
    system = p471.polynomial_system()
    point, residual, jacobian = p471.numerical_root(system)
    radius = Fraction(3, 10**14)
    center = p471.rational_center(point, scale=10**16)
    sines = p471.exact_sine_boxes(scale=10**30)
    preconditioner = p471.rational_matrix(np.linalg.inv(jacobian), scale=10**14)
    normalizer = p471.rational_matrix_from_center(system, center, "normalizer")
    x3 = p471.rational_matrix_from_center(system, center, "x3")
    certificate = {
        "format": "FIN-P480-O180-O181-v1",
        "scope": "declared reduced three-slot channel q=4/5 theta=pi/8",
        "variable_count": 13,
        "parameter_count": 3,
        "equation_count": 13,
        "jacobian_shape": [13, 13],
        "center_scale": 10**16,
        "preconditioner_scale": 10**14,
        "center": [str(value) for value in center],
        "radius": str(radius),
        "sine_boxes": [serialize_interval(value) for value in sines],
        "equations": [serialize_polynomial(value) for value in system["equation_polys"]],
        "jacobian": [
            [serialize_polynomial(value) for value in row]
            for row in system["jacobian_polys"]
        ],
        "preconditioner": serialize_matrix(preconditioner),
        "normalizer_center": serialize_matrix(normalizer),
        "X3_center": serialize_matrix(x3),
        "sylvester_test_margin": "1/100",
        "normalizer_entry_lipschitz": 5,
        "X3_entry_lipschitz": 1,
        "matrix_dimension": 8,
        "lambda_coordinate": 3,
        "floating_locator_residual_infinity_norm": float(np.linalg.norm(residual, np.inf)),
        "boundaries": {
            "network_used": False,
            "laboratory_data_used": False,
            "external_audit_used": False,
            "selector_QW_2191": "open",
            "legacy_strict_bridge": "incomplete_no_role_transfer",
            "physical_status": "dimensionless mathematical operator certificate",
        },
    }
    CERTIFICATE.write_text(
        json.dumps(certificate, indent=2, ensure_ascii=False)+"\n",
        encoding="utf-8",
    )
    print(CERTIFICATE)


if __name__ == "__main__":
    main()
