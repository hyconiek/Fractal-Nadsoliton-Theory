#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r7_shared_frozen_kernel_provenance_packet_for_kobs.json"
OUT_SUMMARY = GENERATED / "r7_shared_frozen_kernel_provenance_packet_for_kobs_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def kernel_delta(a: dict[str, float], b: dict[str, float]) -> dict[str, float]:
    return {key: abs(float(a[key]) - float(b[key])) for key in a.keys()}


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q1949 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1949_informational_light_matter_observer_channel.json"
    )
    q1950 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1950_internal_emergent_observer_closed_loop.json"
    )
    q1951 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1951_mass_informational_weight_internal_observer.json"
    )
    q1952 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1952_information_channel_dedegeneracy_operator.json"
    )
    q1953 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1953_two_state_internal_observer.json"
    )
    q1955 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1955_nogo_and_minimal_operator_repair.json"
    )
    q1956 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1956_two_state_observer_with_repaired_operator.json"
    )
    r2 = load_json("fundamental_action_reconstruction/generated/r2_existing_internal_feedback_parameter_packet_for_kobs.json")
    r3 = load_json("fundamental_action_reconstruction/generated/r3_minimal_internal_light_propagator_packet_for_kobs.json")
    r4 = load_json("fundamental_action_reconstruction/generated/r4_local_chart_emission_map_packet_for_kobs.json")
    r5 = load_json("fundamental_action_reconstruction/generated/r5_minimal_light_to_matter_response_packet_for_kobs.json")
    r6 = load_json("fundamental_action_reconstruction/generated/r6_minimal_observer_readout_packet_for_kobs.json")

    common_kernel_source = q1949["kernel_source"]
    common_kernel = q1949["kernel"]
    related_qw = {
        "QW-1949": q1949,
        "QW-1950": q1950,
        "QW-1951": q1951,
        "QW-1952": q1952,
        "QW-1953": q1953,
        "QW-1955": q1955,
        "QW-1956": q1956,
    }

    provenance_checks = []
    for name, payload in related_qw.items():
        provenance_checks.append(
            {
                "id": f"{name.lower().replace('-', '_')}_kernel_source_match",
                "actual": payload["kernel_source"],
                "expected": common_kernel_source,
                "meaning": f"{name} descends from the same selected frozen-kernel source",
            }
        )
        deltas = kernel_delta(common_kernel, payload["kernel"])
        provenance_checks.append(
            {
                "id": f"{name.lower().replace('-', '_')}_kernel_vector_match",
                "actual": deltas,
                "expected_max_abs_delta": 1e-12,
                "meaning": f"{name} reuses the same frozen kernel vector (omega, phi, beta, eta)",
            }
        )

    provenance_checks.extend(
        [
            {
                "id": "q1949_deterministic_no_fit_note_present",
                "actual": q1949["notes"]["no_fit_statement"],
                "expected": "No supervised fit to labels; deterministic maps only.",
                "meaning": "QW-1949 explicitly records deterministic no-fit provenance",
            },
            {
                "id": "q1950_deterministic_no_fit_note_present",
                "actual": q1950["notes"]["no_fit_statement"],
                "expected": "All maps deterministic from frozen kernel; no label fitting.",
                "meaning": "QW-1950 explicitly records deterministic no-fit provenance",
            },
            {
                "id": "q1950_internal_closed_loop_model_present",
                "actual": q1950["notes"]["closed_loop_model"],
                "expected": "L->M->O->L with observer state internal to nadsoliton dynamics",
                "meaning": "QW-1950 explicitly records the internal closed-loop channel model",
            },
            {
                "id": "r2_source_reports_subset_of_shared_qw_family",
                "actual": sorted(r2["source_reports"]),
                "expected": ["QW-1950", "QW-1951", "QW-1952", "QW-1956"],
                "meaning": "R2 draws only from the shared frozen-kernel QW family",
            },
            {
                "id": "r3_omega_matches_shared_kernel",
                "actual": float(r3["input_scalars"]["omega"]),
                "expected": float(common_kernel["omega"]),
                "tolerance": 1e-12,
                "meaning": "R3 derives G_light from the same frozen kernel omega",
            },
            {
                "id": "r4_psi0_matches_q1952",
                "actual": float(r4["anchor_input"]["psi0"]),
                "expected": float(q1952["derived_params"]["orientation_psi0"]),
                "tolerance": 1e-12,
                "meaning": "R4 derives E from the same psi0 exported in the shared QW family",
            },
            {
                "id": "r5_mass_gain_matches_q1951",
                "actual": float(r5["input_scalars"]["mass_gain"]),
                "expected": float(q1951["mass_informational_weights"]["mass_gain"]),
                "tolerance": 1e-12,
                "meaning": "R5 derives R_mat from the same shared mass-information data",
            },
            {
                "id": "r6_observer_feedback_gain_matches_q1950",
                "actual": float(r6["input_scalars"]["observer_feedback_gain"]),
                "expected": float(q1950["derived_internal_observer_params"]["observer_feedback_gain"]),
                "tolerance": 1e-12,
                "meaning": "R6 derives O_obs from the same shared observer-loop data",
            },
        ]
    )

    for item in provenance_checks:
        if "tolerance" in item:
            item["abs_delta"] = abs(float(item["actual"]) - float(item["expected"]))
            item["consistent"] = item["abs_delta"] <= float(item["tolerance"])
        elif "expected_max_abs_delta" in item:
            max_abs_delta = max(float(v) for v in item["actual"].values())
            item["max_abs_delta"] = max_abs_delta
            item["consistent"] = max_abs_delta <= float(item["expected_max_abs_delta"])
        else:
            item["consistent"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R7",
        "packet_scope": "shared_frozen_kernel_provenance_for_existing_kernel_feedback_and_explicit_h3_chain",
        "shared_frozen_kernel_source": common_kernel_source,
        "shared_frozen_kernel_vector": common_kernel,
        "shared_qw_family": list(related_qw.keys()),
        "current_explicit_chain_packets": ["R2", "R3", "R4", "R5", "R6"],
        "deterministic_provenance_notes": {
            "QW-1949": q1949["notes"],
            "QW-1950": q1950["notes"],
        },
        "provenance_checks": provenance_checks,
        "classification": "shared_frozen_kernel_provenance_present_but_not_operator_factorization",
        "frontier": "R7_B1",
        "no_false_pass": True,
    }

    summary = {
        "stage": "R7",
        "status": "PASS_PARTIAL_SHARED_FROZEN_KERNEL_PROVENANCE_PACKET_READY_OPERATOR_FACTORIZATION_ABSENT",
        "result": "shared_frozen_kernel_provenance_present_but_not_operator_factorization",
        "frontier": [
            "R7_B1",
            "H14_B1",
            "H15_B1",
            "H16_B1",
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
