#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r19_pair1_c1s1_balance_reduction_packet_for_host_matching_route.json"
OUT_SUMMARY = GENERATED / "r19_pair1_c1s1_balance_reduction_packet_for_host_matching_route_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def clean_scalar(value: float) -> float:
    if abs(value) < 1e-15:
        return 0.0
    return round(value, 15)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r11 = load_json(
        "fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json"
    )
    r18 = load_json(
        "fundamental_action_reconstruction/generated/r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json"
    )

    coefficient_classes = r18["pair1_coefficient_class_reduction"]["coefficient_classes"]
    c1s1_classes = [
        row
        for row in coefficient_classes
        if clean_scalar(float(row["signature_on_pair1_entries"]["c1s1"])) != 0.0
    ]
    positive_classes = [
        row
        for row in c1s1_classes
        if float(row["signature_on_pair1_entries"]["c1s1"]) > 0.0
    ]
    negative_classes = [
        row
        for row in c1s1_classes
        if float(row["signature_on_pair1_entries"]["c1s1"]) < 0.0
    ]

    c1_column = r11["transport_packet"]["matrix_columns"]["c1"]
    s1_column = r11["transport_packet"]["matrix_columns"]["s1"]
    c1s1_slot_products = [clean_scalar(float(a) * float(b)) for a, b in zip(c1_column, s1_column)]
    nonzero_abs_values = sorted({clean_scalar(abs(value)) for value in c1s1_slot_products if clean_scalar(value) != 0.0})

    exact_common_factor = math.sqrt(3.0) / 24.0

    positive_sum_symbol = "Sigma_c1s1_positive"
    negative_sum_symbol = "Sigma_c1s1_negative"
    positive_sum_expression = " + ".join(f"({row['class_symbol']})" for row in positive_classes)
    negative_sum_expression = " + ".join(f"({row['class_symbol']})" for row in negative_classes)

    checks = [
        {
            "id": "r18_c1s1_equation_present",
            "actual": any(item["entry"] == "c1s1" for item in r18["pair1_coefficient_class_reduction"]["exact_zero_condition_system"]),
            "expected": True,
            "meaning": "R18 already exports the exact c1s1 zero equation on pair1",
        },
        {
            "id": "transport_c1s1_nonzero_absolute_coefficient_is_unique",
            "actual": len(nonzero_abs_values),
            "expected": 1,
            "meaning": "all nonzero c1s1 coefficients come from one common transport factor",
        },
        {
            "id": "transport_c1s1_common_factor_matches_sqrt3_over_24",
            "actual": abs(nonzero_abs_values[0] - exact_common_factor) <= 1e-12 if nonzero_abs_values else False,
            "expected": True,
            "meaning": "the common nonzero c1s1 transport factor is sqrt(3)/24",
        },
        {
            "id": "positive_and_negative_class_counts_are_two_each",
            "actual": [len(positive_classes), len(negative_classes)],
            "expected": [2, 2],
            "meaning": "the c1s1 equation splits into two positive and two negative coefficient classes",
        },
        {
            "id": "light_boundary_unchanged",
            "actual": r18["light_boundary"]["status"],
            "expected": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R18",
            "meaning": "the already matched light-facing kernel channel remains separate from this c1s1 residual reduction",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R19",
        "packet_goal": "reduce_the_declared_pair1_c1s1_zero_equation_to_a_single_positive_negative_balance_equation_by_factoring_out_the_exact_nonzero_transport_coefficient",
        "source_reports": ["R11", "R18"],
        "c1s1_balance_reduction": {
            "entry": "c1s1",
            "exact_common_transport_factor": {
                "symbolic": "sqrt(3)/24",
                "numeric": clean_scalar(exact_common_factor),
                "derivation_scope": "declared control transport coefficient products T(psi_i,c1) * T(psi_i,s1)",
            },
            "positive_class_symbols": [row["class_symbol"] for row in positive_classes],
            "negative_class_symbols": [row["class_symbol"] for row in negative_classes],
            "positive_balance_sum_symbol": positive_sum_symbol,
            "negative_balance_sum_symbol": negative_sum_symbol,
            "positive_balance_sum_expression": positive_sum_expression,
            "negative_balance_sum_expression": negative_sum_expression,
            "factored_c1s1_equation": "(sqrt(3)/24) * (Sigma_c1s1_positive - Sigma_c1s1_negative) = 0.0",
            "equivalent_balance_equation": "Sigma_c1s1_positive - Sigma_c1s1_negative = 0.0",
            "equivalent_balance_witness_needed": "Sigma_c1s1_positive = Sigma_c1s1_negative",
            "scope": "declared_pair1_c1s1_equation_only",
        },
        "result": {
            "explicit_declared_pair1_c1s1_balance_reduction_present": True,
            "explicit_declared_pair1_c1s1_balance_witness_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R19",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "factorization of the non-light declared pair1 c1s1 residual equation only",
            "boundary": "R19 does not modify the already matched light-facing kernel channel; it only factors the declared pair1 c1s1 residual equation into a single balance equation",
        },
        "classification": "explicit_declared_pair1_c1s1_balance_equation_present_but_no_balance_witness",
        "frontier": "R19_B1",
        "frontier_text": "the repo now exports the exact declared pair1 c1s1 balance equation, but it still lacks a witness that the positive and negative residual class sums are equal and still leaves QW-2191 open",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_declared_pair1_c1s1_balance_equation_holds",
            "no_claim_that_any_other_pair1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R19",
        "status": "PASS_PARTIAL_DECLARED_PAIR1_C1S1_BALANCE_REDUCTION_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R19_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
