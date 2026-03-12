#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_R18 = (
    GENERATED
    / "r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json"
)

OUT_JSON = GENERATED / "p426_current_strict_canonical_local_diagonal_mode2_selector_defect_audit_probe.json"
OUT_SUMMARY = (
    GENERATED / "p426_current_strict_canonical_local_diagonal_mode2_selector_defect_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def linear_combo(coeffs: list[str], terms: list[str]) -> str:
    assert len(coeffs) == len(terms)
    pieces: list[tuple[str, str]] = []
    for coeff, term in zip(coeffs, terms, strict=True):
        if coeff == "0":
            continue
        sign = "+"
        abs_coeff = coeff
        if coeff.startswith("-"):
            sign = "-"
            abs_coeff = coeff[1:]
        if abs_coeff == "1":
            expr = f"({term})"
        else:
            expr = f"({abs_coeff})*({term})"
        pieces.append((sign, expr))
    if not pieces:
        return "0.0"
    out = pieces[0][1]
    if pieces[0][0] == "-":
        out = f"-{out}"
    for sign, expr in pieces[1:]:
        out += f" {sign} {expr}"
    return out


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r18 = load_json(IN_R18)
    classes = {row["class_symbol"]: row for row in r18["pair1_coefficient_class_reduction"]["coefficient_classes"]}

    # Order corresponds to k=0..5 pairs (psi_k, psi_{k+6}).
    class_order = [
        "Sigma_psi0_psi6",
        "Sigma_psi1_psi7",
        "Sigma_psi2_psi8",
        "Sigma_psi3_psi9",
        "Sigma_psi4_psi10",
        "Sigma_psi5_psi11",
    ]

    S_expr: list[str] = []
    S_rows: list[dict[str, Any]] = []
    for symbol in class_order:
        row = classes[symbol]
        S_expr.append(row["residual_sum_expression"])
        S_rows.append(
            {
                "class_symbol": symbol,
                "carrier_slots": row["carrier_slots"],
                "residual_sum_expression": row["residual_sum_expression"],
            }
        )

    # For n=12, mode-2 phase factors are exp(i*pi*k/3) for k=0..5.
    # Re coefficients: [1, 1/2, -1/2, -1, -1/2, 1/2]
    # Im coefficients: [0, sqrt(3)/2, sqrt(3)/2, 0, -sqrt(3)/2, -sqrt(3)/2]
    re_coeffs = ["1", "1/2", "-1/2", "-1", "-1/2", "1/2"]
    im_coeffs = ["0", "sqrt(3)/2", "sqrt(3)/2", "0", "-sqrt(3)/2", "-sqrt(3)/2"]

    re_f2 = linear_combo(re_coeffs, S_expr)
    im_f2 = linear_combo(im_coeffs, S_expr)

    # Pair1 diagonal O(2)-cut signature (a1-d1, b1) from N466:
    # a1-d1 = (2/n) Re(F2), b1 = (1/n) Im(F2), here n=12.
    a1_minus_d1 = f"(1/6)*({re_f2})"
    b1 = f"(1/12)*({im_f2})"

    artifact = {
        "stage": "P426",
        "goal": "materialize_the_exact_mode2_F2_defect_of_the_exported_canonical_local_diagonal_residual_sector_on_n12_as_a_checkable_linear_combination_of_the_six_R18_opposite_pair_sums",
        "inputs": {
            "r18_pair1_residual_declared_pullback_coefficient_class_reduction": str(IN_R18.relative_to(REPO)),
        },
        "mode2_defect": {
            "n": 12,
            "pair_sums_Sk": S_rows,
            "F2_definitions": {
                "F2": "F2(d) := sum_{i=0..11} d_i * exp(i*4*pi*i/12)",
                "reduction": "F2(d) = sum_{k=0..5} S_k * exp(i*pi*k/3) where S_k := d_k + d_{k+6}",
                "Re(F2)_coeffs_on_Sk": re_coeffs,
                "Im(F2)_coeffs_on_Sk": im_coeffs,
            },
            "F2_in_terms_of_R18_pair_sums": {
                "Re_F2_expression": re_f2,
                "Im_F2_expression": im_f2,
            },
            "pair1_o2_cut_signature_via_N466": {
                "Delta1_components": ["a1_minus_d1", "b1"],
                "a1_minus_d1_expression": a1_minus_d1,
                "b1_expression": b1,
                "note": "Diagonal sector breaks O(2) on pair1 iff (a1_minus_d1, b1) != (0,0), equivalently iff F2(d) != 0.",
            },
        },
        "result": {
            "explicit_mode2_defect_expression_present": True,
            "evaluated_numeric_mode2_defect_present": False,
            "strict_core_pair1_o2_cut_witness_present": False,
        },
        "frontier": "P426_B1",
        "frontier_text": "explicit mode-2 defect expression is now exported as a six-class linear combination, but no strict-derived coefficient values exist to decide nonzero; therefore no strict-core O(2)-cut witness is obtained",
        "hard_limits": [
            "no_claim_that_F2_is_nonzero",
            "no_claim_that_the_canonical_local_diagonal_sector_breaks_O2_on_pair1",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_of_strict_core_theta_export",
            "no_claim_of_ToE_closure",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P426",
        "status": "PASS_PARTIAL_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_EXPORTED_VALUE_OPEN",
        "result": "explicit_mode2_defect_expression_present_but_not_evaluable_on_current_repo_exports",
        "frontier": ["P426_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
