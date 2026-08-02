#!/usr/bin/env python3
"""FIN Program P486: exact orientation-premise certificate on the P473 box."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
from typing import Any

import sympy as sp

import fin_phase_exact_algebra as algebra


ROOT = Path(__file__).resolve().parent
RESULT_PATH = ROOT / "FIN_Program_486_Results.json"


def interval_record(value: algebra.Interval) -> dict[str, Any]:
    return {
        "lower": str(value[0]),
        "upper": str(value[1]),
        "lower_float": float(value[0]),
        "upper_float": float(value[1]),
        "strictly_positive": value[0] > 0,
        "strictly_negative": value[1] < 0,
        "excludes_zero": value[0] > 0 or value[1] < 0,
    }


def main() -> None:
    context = algebra.exact_context()
    pivot = algebra.best_active_pivot(context)
    tangent = algebra.tangent_consistency(context)
    box = algebra.certificate_box(context)
    all_symbols = context["system"]["variables"] + (context["alpha"],)

    determinant_c = sp.expand(context["C"].det())
    determinant_c_interval = algebra.polynomial_interval(
        determinant_c, all_symbols, box
    )
    dual_b = context["system"]["variables"][5]
    L_symbol = context["system"]["variables"][3]
    branch_pivot_determinant = sp.expand(pivot["determinant"].subs({
        dual_b: L_symbol-sp.Rational(36, 125)*context["alpha"]
    }))
    pivot_interval = algebra.polynomial_interval(
        branch_pivot_determinant, all_symbols, box
    )
    reference = tangent["reference"]
    reference_coefficient = tangent["coefficients"][reference]
    reference_interval = algebra.polynomial_interval(
        reference_coefficient, all_symbols, box
    )

    A, B, u, L, a, b, c, d, e, f, g, h, i = context["system"]["variables"]
    l_minus_b: algebra.Interval = (
        box[L][0]-box[b][1], box[L][1]-box[b][0]
    )
    normalizer_b = box[B]
    alpha_interval = box[context["alpha"]]

    # Recompute the exact one-dimensional standard-sector factorization.
    system = context["system"]
    eye = sp.eye(8)
    antisymmetric_plus = (
        eye[:, 1]-eye[:, 2]+eye[:, 6]-eye[:, 5]
    )/2
    residual = sp.expand(
        system["x3"]*system["normalizer"]*system["x3"]
        - system["delta"]*system["normalizer"]*system["delta"]/4
    )
    standard_residual = sp.factor(
        (antisymmetric_plus.T*residual*antisymmetric_plus)[0]
    )
    expected = -B*(-125*L+125*b+36*system["sines"][0]) * (
        125*L-125*b+36*system["sines"][0]
    )/15625
    standard_factorization_exact = sp.expand(standard_residual-expected) == 0

    orientation_paid = (
        determinant_c_interval[0] > 0
        and (pivot_interval[0] > 0 or pivot_interval[1] < 0)
        and (reference_interval[0] > 0 or reference_interval[1] < 0)
        and l_minus_b[0] > 0
        and normalizer_b[0] > 0
        and alpha_interval[0] > 0
        and standard_factorization_exact
    )
    result = {
        "metadata": {
            "program": "P486",
            "execution_mode": "local exact rational interval and symbolic audit",
            "network_used": False,
            "laboratory_data_used": False,
            "external_audit_used": False,
        },
        "status": (
            "[Computer-assisted proof] all determinant-sign, positive-branch, "
            "pivot-localization, and rho-reference premises required by P485/P487 "
            "are certified on the exact P473 box"
            if orientation_paid else
            "[Open] at least one P486 orientation/localization premise was not certified"
        ),
        "orientation_premises_paid": orientation_paid,
        "det_C": interval_record(determinant_c_interval),
        "active_pivot_rows": list(pivot["rows"]),
        "active_pivot_determinant": interval_record(pivot_interval),
        "rho_reference_tangent_orbit": reference,
        "rho_reference_coefficient": interval_record(reference_interval),
        "normalizer_B_coordinate": interval_record(normalizer_b),
        "L_minus_b": interval_record(l_minus_b),
        "alpha_sin_pi_over_8": interval_record(alpha_interval),
        "standard_sector_factorization": str(standard_residual),
        "standard_sector_factorization_exact": standard_factorization_exact,
        "positive_branch_theorem": (
            "On the P473 positive box, B>0, L-b>0, and alpha>0. The exact "
            "standard residual therefore selects 125*(L-b)-36*alpha=0 and rejects "
            "the negative factor."
        ),
        "orientation_theorem": (
            "P473 gives X_->0 and N0>0. Since det(C)>0, "
            "Bmap=X_-^{-1} C^T/2 has positive determinant. The exact active "
            "Riccati identity gives Bmap*N0*Bmap^T=N0, hence det(Bmap)^2=1; "
            "therefore det(Bmap)=1 exactly."
        ),
        "necessity": (
            "The pivot interval licenses rational elimination of A/B/u; the "
            "reference-coefficient interval licenses rho=-constant/coefficient; "
            "the orientation proof licenses the odd-dimensional fixed-axis lemma."
        ),
        "boundary": (
            "P486 pays premises only. It does not prove that the fixed axis lies "
            "in the causal equal-endpoint plane; that is P485. It supplies no "
            "selector, physical unit, apparatus, or laboratory evidence."
        ),
        "new_object": "O188 P473 Orientation-and-Localization Premise Certificate",
    }
    RESULT_PATH.write_text(
        json.dumps(algebra.json_ready(result), indent=2, ensure_ascii=False)+"\n",
        encoding="utf-8",
    )
    print(json.dumps(algebra.json_ready(result), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
