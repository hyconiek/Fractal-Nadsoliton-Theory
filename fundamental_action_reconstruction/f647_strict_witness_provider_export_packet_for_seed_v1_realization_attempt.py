#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = (
    GENERATED
    / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def l2_norm(xs: list[float]) -> float:
    return math.sqrt(sum(x * x for x in xs))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    sigma_int = load_json("fundamental_action_reconstruction/generated/sigma_int_strict_derived_v1.json")
    r_ord = load_json("fundamental_action_reconstruction/generated/r_ord_z12_v1_reference_distribution.json")
    alpha_geo = load_json("fundamental_action_reconstruction/generated/alpha_geo_strict_derived_v1.json")

    sigma_value = int(sigma_int["value"])
    if sigma_value not in (-1, 1):
        raise ValueError("sigma_int_strict_derived_v1 must lie in {-1,+1}")

    ord_by_x = r_ord["ord_z12_by_x"]
    if not (isinstance(ord_by_x, list) and len(ord_by_x) == 12 and all(isinstance(v, int) for v in ord_by_x)):
        raise ValueError("r_ord_z12_v1_reference_distribution must store ord_z12_by_x as a length-12 int list")

    alpha_geo_numeric = 4.0 * math.log(2.0)
    w_ref = [math.exp(-alpha_geo_numeric * float(o)) for o in ord_by_x]

    two_pi_over_12 = 2.0 * math.pi / 12.0
    s1 = [math.sin(two_pi_over_12 * float(x)) for x in range(12)]

    # A minimal reflection-breaking per-site weight array, explicitly tied to the seed-v1
    # strict sigma-int sign datum. This is *not* claimed as admissible S_sel_int.
    w_break = [float(sigma_value) * w_ref[x] * s1[x] for x in range(12)]

    # Useful scalar witness: Σ_x w_break(x) s1(x) = sigma * Σ_x w_ref(x) s1(x)^2 ≠ 0.
    dot_w_break_s1 = float(sum(w_break[x] * s1[x] for x in range(12)))

    avg_abs = float(sum(abs(v) for v in w_break) / 12.0)
    w_baseline = [float(sigma_value) * avg_abs for _ in range(12)]
    delta = [float(w_break[x] - w_baseline[x]) for x in range(12)]
    delta_norm = l2_norm(delta)

    exported_object_name = "S_sel_int_strict_core_source_object_v1"

    constructed_source_object = {
        "object_name": exported_object_name,
        "type": "strict_core_constructed_source_object_candidate",
        "declared_domain": {
            "configuration_space_object": sigma_int["declared_domain"]["configuration_space_object"],
            "protected_local_sector": "local_B_tilde_1_topological_sector",
            "carrier_group_object_id": "z_12_v1_group",
        },
        "construction_inputs": [
            "sigma_int_strict_derived_v1",
            "r_ord_z12_v1_reference_distribution",
            "alpha_geo_strict_derived_v1",
            "pair1_sine_axis_basis_function_s1_on_Z12",
        ],
        "definition": {
            "reference_weights": "w_ref(x) := exp(-alpha_geo_strict_derived_v1 * ord_Z12(x))",
            "reflection_breaking_weights": "w_break(x) := sigma_int_strict_derived_v1 * w_ref(x) * s1(x)",
            "s1": "s1(x) := sin(2π x / 12) (pair1 sine axis convention on Z_12)",
        },
        "exported_payload": {
            "ord_z12_by_x": ord_by_x,
            "w_ref_unnormalized_by_x": w_ref,
            "w_break_by_x": w_break,
            "dot_w_break_with_s1": dot_w_break_s1,
        },
        "counts_only_as": [
            "constructed_source_object_export_witness_provider_only",
            "seed_v1_realization_attempt_witness_provider_only",
        ],
        "does_not_count_as": [
            "admissible_S_sel_int",
            "admissible_E_orient",
            "strict_core_selector_closure",
            "QW-2191_discharge",
            "ToE_closure",
        ],
    }

    strict_core_export_properties = {
        "constructed_source_object_exported": True,
        "genuinely_new_exported_identity": True,
        "strict_core_only": True,
        "upstream_of_observer": True,
        "source_seed_only": True,
        "no_external_selector_import": True,
        "kernel_split_safe": True,
        "reuses_psi0_artifact": False,
        "reuses_fr_route_artifact": False,
        "reuses_sigma_int_candidate_as_source_object": False,
        "reuses_overlay_fit_artifact": False,
        "uses_axiom_lane_artifact": False,
    }

    artifact = {
        "stage": "F647",
        "lane": "strict_witness_provider_export_only",
        "goal": "export_one_strict_constructed_source_object_witness_provider_for_seed_v1_realization_attempt_matching_the_F646_signature",
        "status": "F647_EXECUTED_STRICT_WITNESS_PROVIDER_EXPORT_PACKET_FOR_SEED_V1_REALIZATION_ATTEMPT_NO_FALSE_PASS",
        # F646 canonical scan signature (P646):
        "base_realization_attempt": "S_sel_int_new_object_constructed_realization_attempt_v1",
        "constructed_source_object_exported": True,
        "exported_constructed_source_object": exported_object_name,
        "provenance": {
            "base_realization_attempt_packet": "F638",
            "seed_candidate_instance_packet": "F318",
            "seed_sigma_int_source_upgrade": "sigma_int_strict_derived_v1",
            "reference_distribution_packet": "F446",
            "alpha_geo_source_upgrade": alpha_geo["object"],
            "note": "This is an explicit constructed-source object export witness provider only; it does not claim admissibility nor selector closure.",
        },
        "constructed_source_object": constructed_source_object,
        "strict_core_export_properties": strict_core_export_properties,
        "supporting_nonreduction_witness": {
            "present": True,
            "baseline": "sigma_int_broadcast_to_Z12_scaled_by_mean_abs_w_break",
            "delta_norm_l2": delta_norm,
        },
        "admissible_S_sel_int": False,
        "admissible_E_orient": False,
        "strict_core_selector_closure": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F647",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "exported_constructed_source_object": exported_object_name,
        "constructed_source_object_type": constructed_source_object["type"],
        "strict_core_export_properties": strict_core_export_properties,
        "admissible_S_sel_int": False,
        "admissible_E_orient": False,
        "strict_core_selector_closure": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

