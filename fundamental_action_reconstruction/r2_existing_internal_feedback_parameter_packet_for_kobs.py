#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r2_existing_internal_feedback_parameter_packet_for_kobs.json"
OUT_SUMMARY = GENERATED / "r2_existing_internal_feedback_parameter_packet_for_kobs_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q1950 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1950_internal_emergent_observer_closed_loop.json"
    )
    q1951 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1951_mass_informational_weight_internal_observer.json"
    )
    q1952 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1952_information_channel_dedegeneracy_operator.json"
    )
    q1956 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1956_two_state_observer_with_repaired_operator.json"
    )

    artifact = {
        "stage": "R2",
        "parameter_scope": "internal_light_matter_observer_feedback_parameter_packet",
        "source_reports": ["QW-1950", "QW-1951", "QW-1952", "QW-1956"],
        "observer_loop": {
            "observer_tau": q1950["derived_internal_observer_params"]["observer_tau"],
            "observer_feedback_gain": q1950["derived_internal_observer_params"]["observer_feedback_gain"],
            "observer_feedback_theta": q1950["derived_internal_observer_params"]["observer_feedback_theta"],
            "observer_gain_plus": q1950["derived_internal_observer_params"]["observer_gain_plus"],
            "observer_gain_minus": q1950["derived_internal_observer_params"]["observer_gain_minus"],
        },
        "mass_information": {
            "mass_gain": q1951["mass_informational_weights"]["mass_gain"],
            "heavy_weight_sum": q1951["mass_informational_weights"]["heavy_weight_sum"],
            "light_weight_sum": q1951["mass_informational_weights"]["light_weight_sum"],
            "informational_weights": q1951["mass_informational_weights"]["informational_weights"],
        },
        "anisotropy": {
            "anisotropy_strength": q1952["derived_params"]["anisotropy_strength"],
            "retard_phase": q1952["derived_params"]["retard_phase"],
            "orientation_psi0": q1952["derived_params"]["orientation_psi0"],
        },
        "repaired_two_state": {
            "tau_h": q1956["two_state_params"]["tau_h"],
            "tau_l": q1956["two_state_params"]["tau_l"],
            "g_h": q1956["two_state_params"]["g_h"],
            "g_l": q1956["two_state_params"]["g_l"],
            "a2_even_mode": q1956["repaired_operator_params"]["a2_even_mode"],
            "b1_odd_mode": q1956["repaired_operator_params"]["b1_odd_mode"],
            "b3_odd_mode": q1956["repaired_operator_params"]["b3_odd_mode"],
        },
        "classification": "parameter_packet_only_not_operator_factorization",
        "frontier": "R2_B1",
        "no_false_pass": True,
    }

    summary = {
        "stage": "R2",
        "status": "PASS_PARTIAL_EXISTING_INTERNAL_FEEDBACK_PARAMETER_PACKET_READY_OPERATOR_FACTORIZATION_ABSENT",
        "result": "existing_internal_feedback_parameter_packet_present_but_not_yet_an_operator_level_kobs_factorization",
        "frontier": [
            "R2_B1",
            "H14_B1",
            "H15_B1",
            "T12_B1",
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
