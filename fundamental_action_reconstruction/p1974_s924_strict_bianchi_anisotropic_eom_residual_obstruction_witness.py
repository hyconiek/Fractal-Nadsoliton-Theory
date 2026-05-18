#!/usr/bin/env python3
"""P1974 S924 strict Bianchi anisotropic EOM residual obstruction witness.

This is the corrective follow-up to P1973.  P1973 proved only a scalar local
finite-part transport residual.  Here we test the next honest requirement: do
FRW and diagonal Bianchi-I metric EOM components coincide under that scalar
transport?  The answer is no for a generic anisotropic branch.

The result is therefore an obstruction witness, not a closure proof.
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
OUT = GEN / "p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1973 = load("p1973_s923_strict_frw_bianchi_finite_part_transport_matrix_witness.json")

    H, Hd = sp.symbols("H Hdot", real=True)
    s1, s2, ds1, ds2 = sp.symbols("sigma1 sigma2 dsigma1 dsigma2", real=True)
    s3 = sp.simplify(-s1 - s2)
    ds3 = sp.simplify(-ds1 - ds2)

    h = [H + s1, H + s2, H + s3]
    dh = [Hd + ds1, Hd + ds2, Hd + ds3]

    # Diagonal Bianchi-I Einstein tensor components in an orthonormalized
    # component convention: G00 = H1H2+H1H3+H2H3 and
    # Gii/a_i^2 = -(dot H_j + dot H_k + H_j^2 + H_k^2 + H_j H_k).
    g00_bi = sp.expand(h[0] * h[1] + h[0] * h[2] + h[1] * h[2])
    g00_frw = 3 * H**2
    r00 = sp.factor(sp.simplify(g00_bi - g00_frw))

    spatial_residuals = []
    for i in range(3):
        j, k = [idx for idx in range(3) if idx != i]
        gii_bi = -sp.expand(dh[j] + dh[k] + h[j] ** 2 + h[k] ** 2 + h[j] * h[k])
        gii_frw = -(2 * Hd + 3 * H**2)
        spatial_residuals.append(sp.factor(sp.simplify(gii_bi - gii_frw)))

    residual_vector = [r00] + spatial_residuals
    isotropic_subs = {s1: 0, s2: 0, ds1: 0, ds2: 0}
    isotropic_residual_vector = [sp.simplify(expr.subs(isotropic_subs)) for expr in residual_vector]

    polynomial_nonzero_flags = [bool(expr != 0) for expr in residual_vector]

    sample_points = [
        {H: sp.Rational(1, 1), Hd: sp.Rational(1, 10), s1: sp.Rational(1, 10), s2: sp.Rational(-1, 20), ds1: sp.Rational(1, 100), ds2: sp.Rational(-1, 200)},
        {H: sp.Rational(3, 2), Hd: sp.Rational(-1, 20), s1: sp.Rational(1, 5), s2: sp.Rational(1, 10), ds1: sp.Rational(-1, 50), ds2: sp.Rational(3, 200)},
        {H: sp.Rational(4, 5), Hd: sp.Rational(1, 25), s1: sp.Rational(-1, 8), s2: sp.Rational(1, 16), ds1: sp.Rational(1, 80), ds2: sp.Rational(1, 160)},
    ]

    numeric_rows: list[dict[str, Any]] = []
    for idx, point in enumerate(sample_points):
        exact_values = [sp.simplify(expr.subs(point)) for expr in residual_vector]
        float_values = np.array([float(sp.N(v, 30)) for v in exact_values], dtype=float)
        numeric_rows.append(
            {
                "sample_id": f"anisotropic_eom_sample_{idx}",
                "input": {str(k): str(v) for k, v in point.items()},
                "exact_residual_vector": [str(v) for v in exact_values],
                "scipy_l2_norm": float(la.norm(float_values, ord=2)),
                "nonzero": bool(la.norm(float_values, ord=2) > 1e-12),
            }
        )

    obstruction_detected = any(polynomial_nonzero_flags) and all(row["nonzero"] for row in numeric_rows)
    isotropic_limit_zero = all(expr == 0 for expr in isotropic_residual_vector)

    out = {
        "ledger_id": "P1974_S924_STRICT_BIANCHI_ANISOTROPIC_EOM_RESIDUAL_OBSTRUCTION_WITNESS",
        "packet_id": "P1974",
        "stage_id": "S924",
        "produced_by": "p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1973_present": p1973.get("local_result_kind") == "PASS_ZERO_LOCAL_TRANSPORT_RESIDUAL",
        },
        "background_pair": ["FRW", "diagonal_Bianchi-I"],
        "constraint": "sigma1+sigma2+sigma3=0 and dsigma1+dsigma2+dsigma3=0",
        "component_basis": ["G00", "G11/a1^2", "G22/a2^2", "G33/a3^2"],
        "symbolic_setup": {
            "H1": str(h[0]),
            "H2": str(h[1]),
            "H3": str(h[2]),
            "sigma3": str(s3),
            "dsigma3": str(ds3),
            "G00_BianchiI": str(g00_bi),
            "G00_FRW": str(g00_frw),
            "Gii_over_ai2_rule": "-(dot H_j + dot H_k + H_j^2 + H_k^2 + H_j*H_k)",
        },
        "anisotropic_eom_residual_vector": [str(expr) for expr in residual_vector],
        "polynomial_nonzero_flags": polynomial_nonzero_flags,
        "isotropic_limit_residual_vector": [str(expr) for expr in isotropic_residual_vector],
        "isotropic_limit_zero": isotropic_limit_zero,
        "numeric_probe_table": numeric_rows,
        "gatekeeper_checks": {
            "p1973_scalar_transport_available": p1973.get("local_result_kind") == "PASS_ZERO_LOCAL_TRANSPORT_RESIDUAL",
            "generic_anisotropic_residual_nonzero": obstruction_detected,
            "frw_isotropic_limit_zero": isotropic_limit_zero,
            "numeric_samples_nonzero": all(row["nonzero"] for row in numeric_rows),
        },
        "result_kind": "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "SCALAR_TRANSPORT_INSUFFICIENT_FOR_TENSOR_EOM",
            "NONZERO_BIANCHI_I_ANISOTROPIC_EINSTEIN_RESIDUAL",
            "FULL_VARIATIONAL_BACKGROUND_INDEPENDENCE_NOT_CLOSED",
        ],
        "theorem_export": {
            "negative_statement": "The P1973 scalar finite-part transport does not imply full FRW/Bianchi-I metric EOM covariance.  For a generic diagonal Bianchi-I anisotropy with sum_i sigma_i=0, the component residual vector G_BI-G_FRW is polynomially nonzero, while it vanishes in the isotropic FRW limit.",
            "required_new_provider": "A strict anisotropic stress/source tensor or a genuinely tensorial transport connection must cancel/map the displayed residual components before background-independence can be promoted.",
            "not_licensed": [
                "global background-independence closure",
                "PO2/DELTA_BG_Yf sufficiency",
                "PO3 nonempty global admissible class",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1974 is a negative obstruction witness: the scalar transport pass from P1973 must not be upgraded to tensor-level background-independence without an additional anisotropic tensor provider.",
        "next_honest_step": "Construct or rule out a strict anisotropic stress/source tensor Pi_ij^strict whose FRW limit is zero and whose Bianchi-I components cancel the exported residual vector componentwise.",
        "lay_explanation": "Poprzedni most działał dla jednego wspólnego mnożnika.  Teraz sprawdziliśmy pełniejsze równania geometrii i wyszło, że anizotropia Bianchi-I dodaje realne dodatkowe składniki.  To znaczy, że teoria potrzebuje jeszcze jawnego tensorowego źródła albo mocniejszego transportu, zanim można mówić o pełnej niezależności od tła.",
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
