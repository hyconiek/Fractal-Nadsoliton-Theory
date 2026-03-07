#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "f2_strict_gate_kernel_provenance_and_far_input_classification_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f2_strict_gate_kernel_provenance_and_far_input_classification_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def contains_all(text: str, parts: list[str]) -> bool:
    return all(part in text for part in parts)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f1 = load_json(
        "fundamental_action_reconstruction/generated/f1_canonical_informational_fractal_substrate_parameter_packet_summary.json"
    )
    p47 = load_json(
        "fundamental_action_reconstruction/generated/p47_legacy_ontological_kernel_to_strict_gate_kernel_bridge_probe_summary.json"
    )
    n50 = load_json(
        "fundamental_action_reconstruction/generated/n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
    )
    qw2038 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2038_DERIVATION_COMPATIBLE_KERNEL_REFREEZE_SCAN.md"
    )
    qw2039 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2039_DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE.md"
    )
    qw2041 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2041_CANONICAL_REFROZEN_REPARAMETERIZATION_AUDIT.md"
    )
    qw2048 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.md"
    )
    qw2049 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.md"
    )
    qw2064 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2064_MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE.md"
    )

    checks_spec = [
        {
            "id": "f1_legacy_parameter_layer_restored",
            "actual": f1["status"],
            "expected": "F1_EXECUTED_CANONICAL_INFORMATIONAL_FRACTAL_SUBSTRATE_PARAMETER_PACKET_NO_FALSE_PASS",
            "meaning": "F1 restores the canonical ontological parameter layer D_f/alpha_geo/beta_tors inside FAR",
        },
        {
            "id": "p47_no_rigorous_bridge_present",
            "actual": p47["status"],
            "expected": "CURRENT_REPO_EXPORTS_LEGACY_AND_STRICT_KERNELS_BUT_NO_RIGOROUS_LEGACY_TO_STRICT_KERNEL_BRIDGE_AFTER_P47",
            "meaning": "P47 confirms that no rigorous legacy-to-strict kernel bridge exists on the current repo state",
        },
        {
            "id": "n50_nonidentification_theorem_discharged",
            "actual": n50["status"],
            "expected": "N50_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_NONIDENTIFICATION_THEOREM_NO_FALSE_PASS",
            "meaning": "N50 turns the bridge absence into a theorem-level current-repo-state nonidentification result",
        },
        {
            "id": "q2038_operational_refreeze_scan_present",
            "actual": contains_all(qw2038, ["DERIVATION_COMPATIBLE_KERNEL_REFREEZE_PASS", "0.185750 / 0.162500 / 1.000000 / 1.800000"]),
            "expected": True,
            "meaning": "QW-2038 selects the later strict tuple through a refreeze scan",
        },
        {
            "id": "q2039_refrozen_baseline_present",
            "actual": contains_all(qw2039, ["DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE_PASS", "0.185750 / 0.162500 / 0.920000 / 1.800000"]),
            "expected": True,
            "meaning": "QW-2039 promotes a refrozen baseline with beta inside derivational CI95",
        },
        {
            "id": "q2041_semantic_drift_explicit",
            "actual": contains_all(qw2041, ["CANONICAL_REFROZEN_REPARAMETERIZATION_FAIL", "DEFINE_EXPLICIT_BRIDGE_OPERATOR_BETWEEN_CANONICAL_AND_EFFECTIVE_SEMANTICS"]),
            "expected": True,
            "meaning": "QW-2041 explicitly confirms semantic drift and asks for a bridge operator",
        },
        {
            "id": "q2048_methodological_phase_repair_present",
            "actual": contains_all(qw2048, ["SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION_PASS", "beta_target: 0.920000", "eta_target: 1.800000"]),
            "expected": True,
            "meaning": "QW-2048 repairs pointwise identifiability around the refrozen target, rather than deriving the old canonical kernel",
        },
        {
            "id": "q2049_working_gate_selection_present",
            "actual": contains_all(qw2049, ["SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS", "0.185750 / 0.162500 / 1.000000 / 1.800000"]),
            "expected": True,
            "meaning": "QW-2049 selects the later strict working gate tuple through micro/Stage-C intersection",
        },
        {
            "id": "q2064_renormalization_interpretation_is_later_layer",
            "actual": contains_all(qw2064, ["MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE_PASS_WITH_WIDE_CI_WARNING", "z_beta_target: 100.000000", "delta_eta_target: 0.800000"]),
            "expected": True,
            "meaning": "QW-2064 is a later renormalization-interpretation layer rather than the original source of the strict kernel tuple",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "F2",
        "lane": "strict_gate_kernel_provenance_and_far_input_classification_current_repo_state_only",
        "goal": "classify_how_the_later_strict_gate_kernel_may_honestly_enter_FAR_after_the_kernel_split_and_nonidentification_results",
        "status": "F2_EXECUTED_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET_NO_FALSE_PASS",
        "reason": "the current repo shows that the later strict gate kernel comes from an operational refreeze and gate-selection chain QW-2038 to QW-2064, while F1 restores the old ontological parameter layer and P47/N50 show that no rigorous bridge currently identifies the strict gate kernel with the old ontological kernel; therefore FAR must classify the strict gate kernel as a later-pipeline operational import or control only, not as a silent ontological replacement for the legacy action-first source layer",
        "strict_gate_kernel_derivation_chain": [
            "QW-2038 derivation-compatible refreeze scan",
            "QW-2039 derivation-compatible refrozen baseline",
            "QW-2041 canonical-vs-refrozen reparameterization fail and bridge obligation",
            "QW-2048 spectral phase-locked pointwise repair",
            "QW-2049 spectral micro/Stage-C intersection selection",
            "QW-2050 freeze bundle",
            "QW-2064 micro-derived renormalization-constants interpretation",
        ],
        "far_input_classification": {
            "legacy_ontological_source_layer": {
                "class": "canonical_ontology_supported_action_first_source_layer",
                "objects": [
                    "K_legacy_ont_or_its_restored_parameter_layer",
                    "D_f",
                    "alpha_geo",
                    "beta_tors",
                ],
            },
            "strict_gate_kernel_layer": {
                "class": "later_pipeline_operational_strict_working_kernel",
                "objects": [
                    "K_strict_gate",
                    "omega=0.18575",
                    "phi=0.16250",
                    "beta=1.0",
                    "eta=1.8",
                ],
            },
            "silent_full_ontological_substitution_disallowed": True,
            "A1_rule": "do_not_silently_replace_the_action_first_ontological_source_layer_with_K_strict_gate",
            "A4_rule": "use_K_strict_gate_only_as_operational_control_or_downstream_consistency_target_unless_a_bridge_is_added",
            "A8_rule": "do_not_silently_inherit_alpha_geo_over_2beta_tors_hierarchy_semantics_into_K_strict_gate",
        },
        "remaining_missing_objects": [
            "explicit_legacy_to_strict_kernel_bridge_operator_or_theorem",
            "explicit_amplitude_absorption_or_normalization_map_for_alpha_geo",
            "explicit_beta_tors_to_beta_eta_translation_rule",
            "explicit_phase_frequency_bridge",
            "explicit_kernel_level_ontological_anchoring_witness_for_K_strict_gate",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F2",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "strict_gate_kernel_derivation_chain": artifact["strict_gate_kernel_derivation_chain"],
        "far_input_classification": artifact["far_input_classification"],
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
