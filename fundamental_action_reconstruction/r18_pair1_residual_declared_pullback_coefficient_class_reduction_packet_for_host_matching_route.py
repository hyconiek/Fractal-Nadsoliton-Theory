#!/usr/bin/env python3
from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json"
)
OUT_SUMMARY = (
    GENERATED
    / "r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def clean_scalar(value: float) -> float:
    if abs(value) < 1e-15:
        return 0.0
    return round(value, 15)


def class_scalar(value: float) -> float:
    if abs(value) < 1e-12:
        return 0.0
    return round(value, 12)


def format_coeff(value: float) -> str:
    cleaned = clean_scalar(value)
    if cleaned == 0.0:
        return "0.0"
    return str(cleaned)


def build_linear_expression(terms: list[tuple[float, str]]) -> str:
    kept_terms: list[str] = []
    for coefficient, symbol in terms:
        cleaned = clean_scalar(coefficient)
        if cleaned == 0.0:
            continue
        if cleaned == 1.0:
            kept_terms.append(symbol)
        elif cleaned == -1.0:
            kept_terms.append(f"-{symbol}")
        else:
            kept_terms.append(f"({cleaned})*({symbol})")
    if not kept_terms:
        return "0.0"

    expression = kept_terms[0]
    for term in kept_terms[1:]:
        if term.startswith("-"):
            expression += " - " + term[1:]
        else:
            expression += " + " + term
    return expression


def build_sum_expression(terms: list[str]) -> str:
    if not terms:
        return "0.0"
    if len(terms) == 1:
        return terms[0]
    return " + ".join(f"({term})" for term in terms)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r15 = load_json(
        "fundamental_action_reconstruction/generated/r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
    )
    r16 = load_json(
        "fundamental_action_reconstruction/generated/r16_explicit_residual_local_diagonal_declared_control_pullback_packet_for_host_matching_route.json"
    )
    r17 = load_json(
        "fundamental_action_reconstruction/generated/r17_explicit_host_side_residual_diagonal_correction_absence_packet_for_host_matching_route.json"
    )

    control_basis = r16["declared_control_pullback_packet"]["control_basis"]
    carrier_basis = r16["declared_control_pullback_packet"]["carrier_basis"]
    coefficient_tensor = r16["declared_control_pullback_packet"]["coefficient_tensor_by_carrier_slot"]
    residual_entries = [
        entry["residual_local_diagonal_entry"]
        for entry in r15["diagonal_decomposition"]["entrywise_rows"]
    ]

    pair1_basis = control_basis[:2]
    coeff_c1c1 = coefficient_tensor[0][0]
    coeff_c1s1 = coefficient_tensor[0][1]
    coeff_s1s1 = coefficient_tensor[1][1]
    coeff_s1c1 = coefficient_tensor[1][0]

    coefficient_classes: OrderedDict[
        tuple[float, float, float], dict[str, Any]
    ] = OrderedDict()
    for index, basis_label in enumerate(carrier_basis):
        signature = (
            class_scalar(coeff_c1c1[index]),
            class_scalar(coeff_c1s1[index]),
            class_scalar(coeff_s1s1[index]),
        )
        if signature not in coefficient_classes:
            coefficient_classes[signature] = {
                "carrier_slots": [],
                "residual_entries": [],
            }
        coefficient_classes[signature]["carrier_slots"].append(basis_label)
        coefficient_classes[signature]["residual_entries"].append(residual_entries[index])

    class_rows: list[dict[str, Any]] = []
    c1c1_terms: list[tuple[float, str]] = []
    c1s1_terms: list[tuple[float, str]] = []
    s1s1_terms: list[tuple[float, str]] = []
    zero_conditions: list[dict[str, Any]] = []

    for ordinal, (signature, payload) in enumerate(coefficient_classes.items(), start=1):
        class_symbol = "Sigma_" + "_".join(payload["carrier_slots"])
        class_row = {
            "class_id": f"CC{ordinal}",
            "class_symbol": class_symbol,
            "carrier_slots": payload["carrier_slots"],
            "signature_on_pair1_entries": {
                "c1c1": signature[0],
                "c1s1": signature[1],
                "s1s1": signature[2],
            },
            "residual_sum_expression": build_sum_expression(payload["residual_entries"]),
        }
        class_rows.append(class_row)
        c1c1_terms.append((signature[0], class_symbol))
        c1s1_terms.append((signature[1], class_symbol))
        s1s1_terms.append((signature[2], class_symbol))

    independent_entry_equations = [
        {
            "entry": "c1c1",
            "exact_linear_combination": build_linear_expression(c1c1_terms),
        },
        {
            "entry": "c1s1",
            "exact_linear_combination": build_linear_expression(c1s1_terms),
        },
        {
            "entry": "s1s1",
            "exact_linear_combination": build_linear_expression(s1s1_terms),
        },
    ]

    for equation in independent_entry_equations:
        zero_conditions.append(
            {
                "entry": equation["entry"],
                "required_zero_equation": f"{equation['exact_linear_combination']} = 0.0",
            }
        )

    checks = [
        {
            "id": "r16_pair1_declared_block_present",
            "actual": r16["result"]["pair1_declared_residual_diagonal_block_present"],
            "expected": True,
            "meaning": "R16 already exports the declared pair1 residual block",
        },
        {
            "id": "r17_host_side_residual_branch_absent",
            "actual": r17["result"]["host_side_residual_diagonal_correction_present"],
            "expected": False,
            "meaning": "R17 already closes the host-side residual correction branch",
        },
        {
            "id": "pair1_cross_entry_is_symmetric",
            "actual": [clean_scalar(value) for value in coeff_c1s1],
            "expected": [clean_scalar(value) for value in coeff_s1c1],
            "meaning": "the pair1 residual block has only three independent entries because the off-diagonal entry is symmetric",
        },
        {
            "id": "transport_induced_coefficient_class_partition_is_complete",
            "actual": sum(len(row["carrier_slots"]) for row in class_rows),
            "expected": len(carrier_basis),
            "meaning": "the coefficient classes partition the full 12-slot carrier",
        },
        {
            "id": "transport_induced_coefficient_class_count_is_six",
            "actual": len(class_rows),
            "expected": 6,
            "meaning": "the declared pair1 residual block reduces to six exact coefficient classes on the current transport",
        },
        {
            "id": "shared_kernel_light_channel_stays_separate",
            "actual": r16["light_boundary"]["status"],
            "expected": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R15_R16",
            "meaning": "the already matched light-facing kernel channel remains separate from this non-light residual diagonal reduction",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R18",
        "packet_goal": "reduce_the_pair1_declared_pullback_of_the_residual_local_diagonal_sector_to_exact_transport_induced_coefficient_classes_and_a_finite_zero_equation_system_without_claiming_that_any_zero_equation_is_satisfied",
        "source_reports": ["R15", "R16", "R17"],
        "pair1_coefficient_class_reduction": {
            "pair1_basis": pair1_basis,
            "source_object": "declared_pair1_block_of_M_control_residual_diag_declared",
            "class_family_kind": "transport_induced_coefficient_classes_on_the_declared_pair1_residual_block",
            "coefficient_classes": class_rows,
            "independent_entry_equations": independent_entry_equations,
            "exact_zero_condition_system": zero_conditions,
            "reduction_scope": "declared_pair1_block_only",
        },
        "result": {
            "explicit_transport_induced_pair1_coefficient_class_reduction_present": True,
            "explicit_zero_condition_system_for_declared_pair1_residual_block_present": True,
            "explicit_zero_witness_for_declared_pair1_residual_block_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R18",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "coefficient-class reduction of the non-light residual local diagonal complement on the declared pair1 block only",
            "boundary": "R18 does not modify the already matched light-facing kernel channel; it only rewrites the residual diagonal pair1 block into a finite coefficient-class zero system",
        },
        "classification": "explicit_declared_pair1_residual_zero_system_present_but_no_zero_witness",
        "frontier": "R18_B1",
        "frontier_text": "the repo now exports the exact finite zero-equation system for the declared pair1 residual block, but it still lacks witnesses that these equations hold and still leaves QW-2191 open",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_any_pair1_zero_equation_is_satisfied",
            "no_claim_that_the_full_residual_declared_pullback_vanishes",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R18",
        "status": "PASS_PARTIAL_PAIR1_RESIDUAL_DECLARED_PULLBACK_COEFFICIENT_CLASS_REDUCTION_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R18_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
