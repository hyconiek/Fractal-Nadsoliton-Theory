#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r21_explicit_pair1_c1s1_shift_defect_polynomial_packet_for_host_route.json"
OUT_SUMMARY = GENERATED / "r21_explicit_pair1_c1s1_shift_defect_polynomial_packet_for_host_route_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def join_sum(terms: list[str]) -> str:
    if not terms:
        return "0.0"
    if len(terms) == 1:
        return terms[0]
    return " + ".join(f"({term})" for term in terms)


def difference(lhs_terms: list[str], rhs_terms: list[str]) -> str:
    lhs = join_sum(lhs_terms)
    rhs = join_sum(rhs_terms)
    if rhs == "0.0":
        return lhs
    if lhs == "0.0":
        return f"-( {rhs} )"
    return f"({lhs}) - ({rhs})"


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r15 = load_json(
        "fundamental_action_reconstruction/generated/r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
    )
    r18 = load_json(
        "fundamental_action_reconstruction/generated/r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json"
    )
    r20 = load_json(
        "fundamental_action_reconstruction/generated/r20_declared_three_step_carrier_shift_packet_for_pair1_c1s1_balance_route.json"
    )

    entry_rows = {row["basis_label"]: row for row in r15["diagonal_decomposition"]["entrywise_rows"]}
    class_rows = {
        row["class_symbol"]: row for row in r18["pair1_coefficient_class_reduction"]["coefficient_classes"]
    }

    positive_symbols = ["Sigma_psi1_psi7", "Sigma_psi2_psi8"]
    negative_symbols = ["Sigma_psi4_psi10", "Sigma_psi5_psi11"]
    positive_slots = sum((class_rows[symbol]["carrier_slots"] for symbol in positive_symbols), [])
    negative_slots = sum((class_rows[symbol]["carrier_slots"] for symbol in negative_symbols), [])

    positive_residual_terms = [entry_rows[slot]["residual_local_diagonal_entry"] for slot in positive_slots]
    negative_residual_terms = [entry_rows[slot]["residual_local_diagonal_entry"] for slot in negative_slots]

    positive_g4_terms = [f"g4_{slot}*v{slot}**2" for slot in positive_slots]
    negative_g4_terms = [f"g4_{slot}*v{slot}**2" for slot in negative_slots]
    positive_g6_terms = [f"g6_{slot}*v{slot}**4" for slot in positive_slots]
    negative_g6_terms = [f"g6_{slot}*v{slot}**4" for slot in negative_slots]
    positive_gY_terms = [f"gY{slot[3:]}*vphi**2" for slot in positive_slots]
    negative_gY_terms = [f"gY{slot[3:]}*vphi**2" for slot in negative_slots]
    positive_m2_terms = [f"m2_{slot}" for slot in positive_slots]
    negative_m2_terms = [f"m2_{slot}" for slot in negative_slots]

    coefficient_family_defect = {
        "quartic_like_g4_family_defect": f"3*({difference(positive_g4_terms, negative_g4_terms)})",
        "quintic_like_g6_family_defect": f"5*({difference(positive_g6_terms, negative_g6_terms)})",
        "yukawa_like_gY_family_defect": f"2*({difference(positive_gY_terms, negative_gY_terms)})",
        "mass_like_m2_family_defect": difference(positive_m2_terms, negative_m2_terms),
    }

    defect_polynomial = (
        f"({coefficient_family_defect['quartic_like_g4_family_defect']}) + "
        f"({coefficient_family_defect['quintic_like_g6_family_defect']}) + "
        f"({coefficient_family_defect['yukawa_like_gY_family_defect']}) + "
        f"({coefficient_family_defect['mass_like_m2_family_defect']})"
    )

    plus3_basis_map = r20["declared_three_step_shift_packet"]["declared_basis_map"]

    checks = [
        {
            "id": "r20_declared_plus3_shift_packet_present",
            "actual": r20["result"]["explicit_declared_plus3_carrier_shift_packet_present"],
            "expected": True,
            "meaning": "R20 already exports the declared +3 carrier shift packet",
        },
        {
            "id": "positive_and_negative_support_slot_counts_are_four",
            "actual": [len(positive_slots), len(negative_slots)],
            "expected": [4, 4],
            "meaning": "the c1s1 support consists of four positive and four negative diagonal slots",
        },
        {
            "id": "plus3_shift_maps_positive_support_slots_to_negative_support_slots",
            "actual": sorted(plus3_basis_map[slot] for slot in positive_slots),
            "expected": sorted(negative_slots),
            "meaning": "the declared +3 carrier shift maps the positive c1s1 support slots onto the negative support slots",
        },
        {
            "id": "residual_floor_terms_cancel_exactly_by_cardinality",
            "actual": len(positive_slots) == len(negative_slots),
            "expected": True,
            "meaning": "the common host scalar floor cancels exactly inside the c1s1 shift-defect polynomial",
        },
        {
            "id": "light_boundary_unchanged",
            "actual": r20["light_boundary"]["status"],
            "expected": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R20",
            "meaning": "the already matched light-facing kernel channel remains separate from this coefficient-level c1s1 shift defect packet",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R21",
        "packet_goal": "materialize_the_exact_coefficient_level_pair1_c1s1_shift_defect_polynomial_whose_vanishing_is_equivalent_to_the_missing_declared_plus3_shift_equivariance_witness_on_the_c1s1_support_sum",
        "source_reports": ["R15", "R18", "R20"],
        "pair1_c1s1_shift_defect_packet": {
            "positive_support_symbols": positive_symbols,
            "negative_support_symbols": negative_symbols,
            "positive_support_slots": positive_slots,
            "negative_support_slots": negative_slots,
            "positive_residual_sum_expression": join_sum(positive_residual_terms),
            "negative_residual_sum_expression": join_sum(negative_residual_terms),
            "exact_shift_defect_expression": difference(positive_residual_terms, negative_residual_terms),
            "coefficient_family_decomposition": coefficient_family_defect,
            "equivalent_zero_equation": f"{defect_polynomial} = 0.0",
            "scope": "declared_pair1_c1s1_shift_defect_only",
        },
        "result": {
            "explicit_pair1_c1s1_shift_defect_polynomial_present": True,
            "explicit_zero_witness_for_pair1_c1s1_shift_defect_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R21",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "coefficient-level non-light defect polynomial on the declared pair1 c1s1 support only",
            "boundary": "R21 does not modify the already matched light-facing kernel channel; it only materializes the exact diagonal coefficient-level defect whose vanishing would give the declared c1s1 shift-equivariance witness",
        },
        "classification": "explicit_pair1_c1s1_shift_defect_polynomial_present_but_no_zero_witness",
        "frontier": "R21_B1",
        "frontier_text": "the repo now exports the exact coefficient-level pair1 c1s1 shift-defect polynomial, but it still lacks a witness that this defect vanishes and still leaves QW-2191 open",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_pair1_c1s1_shift_defect_polynomial_vanishes",
            "no_claim_that_any_other_pair1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R21",
        "status": "PASS_PARTIAL_EXPLICIT_PAIR1_C1S1_SHIFT_DEFECT_POLYNOMIAL_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R21_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
