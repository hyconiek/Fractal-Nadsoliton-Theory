#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r3_minimal_internal_light_propagator_packet_for_kobs.json"
OUT_SUMMARY = GENERATED / "r3_minimal_internal_light_propagator_packet_for_kobs_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q1952 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1952_information_channel_dedegeneracy_operator.json"
    )
    q1955 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1955_nogo_and_minimal_operator_repair.json"
    )
    q1956 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1956_two_state_observer_with_repaired_operator.json"
    )

    omega = float(q1952["kernel"]["omega"])
    retard_phase = float(q1952["derived_params"]["retard_phase"])
    anisotropy_strength = float(q1952["derived_params"]["anisotropy_strength"])
    c_speed = 1.0

    l0 = c_speed * retard_phase / omega
    l_plus = l0 * (1.0 + anisotropy_strength)
    l_minus = l0 * (1.0 - anisotropy_strength)
    lambda_plus = math.cos(omega * l_plus / c_speed)
    lambda_minus = math.cos(omega * l_minus / c_speed)

    consistency_checks = [
        {
            "id": "retard_phase_q1952_q1955_match",
            "actual": retard_phase,
            "expected": float(q1955["minimal_repair_params"]["retard_phase"]),
            "tolerance": 1e-12,
        },
        {
            "id": "retard_phase_q1955_q1956_match",
            "actual": float(q1955["minimal_repair_params"]["retard_phase"]),
            "expected": float(q1956["repaired_operator_params"]["retard_phase"]),
            "tolerance": 1e-12,
        },
    ]
    for item in consistency_checks:
        item["abs_delta"] = abs(item["actual"] - item["expected"])
        item["consistent"] = item["abs_delta"] <= item["tolerance"]

    artifact = {
        "stage": "R3",
        "operator_scope": "explicit_internal_light_propagator_packet",
        "carrier": "L_1=span{ell_+,ell_-}",
        "basis_order": ["ell_+", "ell_-"],
        "basis_role": "ordered_light_eigenchannel_basis_only_not_selector_or_pair_basis",
        "source_reports": ["QW-1952", "QW-1955", "QW-1956"],
        "unit_convention": {"speed_of_signal_c": c_speed},
        "operator_formula": "G_light^(1)=diag(cos(omega*L_+/c),cos(omega*L_-/c))",
        "input_scalars": {
            "omega": omega,
            "retard_phase": retard_phase,
            "anisotropy_strength": anisotropy_strength,
        },
        "derived_path_lengths": {
            "L0": l0,
            "L_plus": l_plus,
            "L_minus": l_minus,
            "path_rule": {
                "L0": "c*retard_phase/omega",
                "L_plus": "L0*(1+anisotropy_strength)",
                "L_minus": "L0*(1-anisotropy_strength)",
            },
        },
        "matrix": [
            [lambda_plus, 0.0],
            [0.0, lambda_minus],
        ],
        "eigenvalues": {
            "lambda_plus": lambda_plus,
            "lambda_minus": lambda_minus,
        },
        "self_adjoint": True,
        "uses_orientation_anchor_in_matrix": False,
        "selector_sector_projected_block_present": False,
        "factorization_status": "not_identified_with_existing_kernel_feedback",
        "pair_projection_status": "absent",
        "consistency_checks": consistency_checks,
        "classification": "explicit_internal_light_propagator_present_but_eigenchannel_only_and_not_kernel_factorized",
        "frontier": "R3_B1",
        "no_false_pass": True,
    }

    summary = {
        "stage": "R3",
        "status": "PASS_PARTIAL_EXPLICIT_INTERNAL_LIGHT_PROPAGATOR_PACKET_READY_FACTORIZATION_ABSENT",
        "result": "explicit_internal_light_propagator_packet_present_but_eigenchannel_only_and_not_yet_factorized_from_existing_kernel_feedback",
        "frontier": [
            "R3_B1",
            "H14_B1",
            "H15_B1",
            "H29_B1",
            "C32_B2",
        ],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

