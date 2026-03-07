#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r20_declared_three_step_carrier_shift_packet_for_pair1_c1s1_balance_route.json"
OUT_SUMMARY = GENERATED / "r20_declared_three_step_carrier_shift_packet_for_pair1_c1s1_balance_route_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def clean_scalar(value: float) -> float:
    if abs(value) < 1e-15:
        return 0.0
    return round(value, 15)


def shift_right(values: list[float], amount: int) -> list[float]:
    amount = amount % len(values)
    return values[-amount:] + values[:-amount]


def same_vector(lhs: list[float], rhs: list[float]) -> bool:
    return all(abs(float(a) - float(b)) <= 1e-12 for a, b in zip(lhs, rhs))


def neg_vector(lhs: list[float], rhs: list[float]) -> bool:
    return all(abs(float(a) + float(b)) <= 1e-12 for a, b in zip(lhs, rhs))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r11 = load_json(
        "fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json"
    )
    r18 = load_json(
        "fundamental_action_reconstruction/generated/r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json"
    )
    r19 = load_json(
        "fundamental_action_reconstruction/generated/r19_pair1_c1s1_balance_reduction_packet_for_host_matching_route.json"
    )

    target_basis = r11["transport_packet"]["target_basis"]
    c1_column = [float(value) for value in r11["transport_packet"]["matrix_columns"]["c1"]]
    s1_column = [float(value) for value in r11["transport_packet"]["matrix_columns"]["s1"]]
    c2_column = [float(value) for value in r11["transport_packet"]["matrix_columns"]["c2"]]
    s2_column = [float(value) for value in r11["transport_packet"]["matrix_columns"]["s2"]]

    plus3_basis_map = {basis: target_basis[(index + 3) % len(target_basis)] for index, basis in enumerate(target_basis)}

    shifted_c1 = shift_right(c1_column, 3)
    shifted_s1 = shift_right(s1_column, 3)
    shifted_c2 = shift_right(c2_column, 3)
    shifted_s2 = shift_right(s2_column, 3)

    class_rows = {
        row["class_symbol"]: row for row in r18["pair1_coefficient_class_reduction"]["coefficient_classes"]
    }
    positive_symbols = r19["c1s1_balance_reduction"]["positive_class_symbols"]
    negative_symbols = r19["c1s1_balance_reduction"]["negative_class_symbols"]

    class_shift_map: dict[str, str] = {}
    for symbol in positive_symbols + negative_symbols:
        carrier_slots = class_rows[symbol]["carrier_slots"]
        shifted_slots = sorted(plus3_basis_map[slot] for slot in carrier_slots)
        for target_symbol, target_row in class_rows.items():
            if sorted(target_row["carrier_slots"]) == shifted_slots:
                class_shift_map[symbol] = target_symbol
                break

    checks = [
        {
            "id": "target_basis_is_12_cycle",
            "actual": len(target_basis),
            "expected": 12,
            "meaning": "the declared transport carrier is the 12-slot octave ring",
        },
        {
            "id": "plus3_shift_maps_c1_to_s1",
            "actual": same_vector(shifted_c1, s1_column),
            "expected": True,
            "meaning": "the declared +3 carrier shift sends c1 to s1 on the transported pair1 plane",
        },
        {
            "id": "plus3_shift_maps_s1_to_minus_c1",
            "actual": neg_vector(shifted_s1, c1_column),
            "expected": True,
            "meaning": "the declared +3 carrier shift sends s1 to -c1 on the transported pair1 plane",
        },
        {
            "id": "plus3_shift_maps_c2_to_minus_c2",
            "actual": neg_vector(shifted_c2, c2_column),
            "expected": True,
            "meaning": "the declared +3 carrier shift preserves the second pair up to sign",
        },
        {
            "id": "plus3_shift_maps_s2_to_minus_s2",
            "actual": neg_vector(shifted_s2, s2_column),
            "expected": True,
            "meaning": "the declared +3 carrier shift preserves the second sine mode up to sign",
        },
        {
            "id": "plus3_shift_maps_positive_c1s1_support_to_negative_support",
            "actual": [class_shift_map[symbol] for symbol in positive_symbols],
            "expected": negative_symbols,
            "meaning": "the declared +3 carrier shift maps the positive c1s1 support classes onto the negative support classes",
        },
        {
            "id": "plus3_shift_maps_negative_c1s1_support_back_to_positive_support",
            "actual": [class_shift_map[symbol] for symbol in negative_symbols],
            "expected": positive_symbols,
            "meaning": "the declared +3 carrier shift maps the negative c1s1 support classes back onto the positive support classes",
        },
        {
            "id": "light_boundary_unchanged",
            "actual": r19["light_boundary"]["status"],
            "expected": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R19",
            "meaning": "the already matched light-facing kernel channel remains separate from this declared carrier-shift packet",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R20",
        "packet_goal": "materialize_the_exact_declared_plus3_carrier_shift_that_flips_the_pair1_c1s1_sign_and_maps_the_positive_support_classes_onto_the_negative_support_classes",
        "source_reports": ["R11", "R18", "R19"],
        "declared_three_step_shift_packet": {
            "carrier_basis": target_basis,
            "declared_basis_map": plus3_basis_map,
            "transported_pair1_action": {
                "S_plus3_c1": "s1",
                "S_plus3_s1": "-c1",
                "consequence_for_c1s1": "c1s1 -> -c1s1",
            },
            "support_class_shift_map": class_shift_map,
            "equivalent_c1s1_witness_reduction": {
                "from_balance_equation": "Sigma_c1s1_positive = Sigma_c1s1_negative",
                "to_shift_equivariance_equation": "S_plus3(Sigma_c1s1_positive) = Sigma_c1s1_positive",
                "because": "S_plus3(Sigma_c1s1_positive) = Sigma_c1s1_negative on the declared support",
            },
            "scope": "declared_pair1_c1s1_support_only",
        },
        "result": {
            "explicit_declared_plus3_carrier_shift_packet_present": True,
            "explicit_declared_plus3_shift_equivariance_witness_for_pair1_c1s1_support_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R20",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "declared +3 carrier-shift structure on the non-light pair1 c1s1 residual support only",
            "boundary": "R20 does not modify the already matched light-facing kernel channel; it only materializes the declared carrier shift that flips the pair1 c1s1 sign",
        },
        "classification": "explicit_declared_plus3_shift_packet_present_but_no_c1s1_shift_equivariance_witness",
        "frontier": "R20_B1",
        "frontier_text": "the repo now exports the exact declared +3 carrier shift for the pair1 c1s1 branch, but it still lacks a witness that the residual c1s1 support sum is invariant under that shift and still leaves QW-2191 open",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_declared_plus3_shift_is_physically_canonical",
            "no_claim_that_the_pair1_c1s1_shift_equivariance_equation_holds",
            "no_claim_that_any_other_pair1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R20",
        "status": "PASS_PARTIAL_DECLARED_PLUS3_CARRIER_SHIFT_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R20_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
