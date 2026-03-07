#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "F1": load_json(
            "fundamental_action_reconstruction/generated/f1_canonical_informational_fractal_substrate_parameter_packet_summary.json"
        ),
        "P47": load_json(
            "fundamental_action_reconstruction/generated/p47_legacy_ontological_kernel_to_strict_gate_kernel_bridge_probe_summary.json"
        ),
    }

    checks_spec = [
        {
            "id": "f1_restores_legacy_parameter_layer",
            "actual": sources["F1"]["status"],
            "expected": "F1_EXECUTED_CANONICAL_INFORMATIONAL_FRACTAL_SUBSTRATE_PARAMETER_PACKET_NO_FALSE_PASS",
            "meaning": "F1 restores the canonical legacy parameter layer inside FAR",
        },
        {
            "id": "p47_confirms_kernel_split_without_bridge",
            "actual": sources["P47"]["status"],
            "expected": "CURRENT_REPO_EXPORTS_LEGACY_AND_STRICT_KERNELS_BUT_NO_RIGOROUS_LEGACY_TO_STRICT_KERNEL_BRIDGE_AFTER_P47",
            "meaning": "P47 confirms that the repo exports both kernels but no rigorous bridge",
        },
        {
            "id": "p47_shared_parameter_bridge_absent",
            "actual": sources["P47"]["bridge_state"]["authoritative_shared_parameter_bridge_present"],
            "expected": False,
            "meaning": "no authoritative shared-parameter bridge object is present",
        },
        {
            "id": "p47_beta_translation_absent",
            "actual": sources["P47"]["bridge_state"]["explicit_beta_tors_to_beta_eta_translation_present"],
            "expected": False,
            "meaning": "no explicit beta_tors to beta/eta translation is present",
        },
        {
            "id": "p47_phase_bridge_absent",
            "actual": sources["P47"]["bridge_state"]["explicit_phase_frequency_bridge_present"],
            "expected": False,
            "meaning": "no explicit phase/frequency bridge is present",
        },
        {
            "id": "p47_ontological_anchoring_witness_absent",
            "actual": sources["P47"]["bridge_state"]["explicit_kernel_level_ontological_anchoring_witness_present"],
            "expected": False,
            "meaning": "the strict kernel still lacks a kernel-level ontological anchoring witness",
        },
    ]

    checks: list[dict[str, Any]] = []
    mismatches: list[str] = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "step": "N50",
            "status": "N50_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_TO_STRICT_KERNEL_BRIDGE_STATE",
            "scope": "current_legacy_to_strict_kernel_bridge_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N50",
            "status": "N50_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_NONIDENTIFICATION_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_to_strict_kernel_bridge_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "legacy_kernel_present": True,
                "strict_gate_kernel_present": True,
                "legacy_parameter_layer_restored": True,
                "rigorous_bridge_present": False,
                "rigorous_nonidentification_on_current_repo_state": True,
                "silent_full_ontological_inheritance_from_legacy_to_strict_is_unsupported": True,
                "strict_kernel_operationally_usable_but_ontologically_underanchored_without_bridge": True,
                "full_closure_pass": False,
            },
            "missing_structure_classes": [
                "explicit_amplitude_normalization_or_absorption_map_explaining_the_loss_of_visible_alpha_geo_between_K_legacy_ont_and_K_strict_gate",
                "explicit_renormalized_damping_map_beta_tors_to_beta_eta_or_equivalent_strict_translation_rule",
                "explicit_phase_frequency_bridge_from_pi_over_4_pi_over_6_to_0p18575_0p16250",
                "explicit_kernel_specific_retained_vs_replaced_partition_theorem_for_legacy_vs_strict_kernel_roles",
                "explicit_kernel_level_ontological_anchoring_witness_showing_that_K_strict_gate_is_an_internal_route_of_one_informational_nadsoliton",
            ],
            "hard_limits": [
                "no_proof_that_the_strict_gate_kernel_is_false",
                "no_proof_that_the_legacy_kernel_is_false",
                "no_proof_that_no_bridge_can_ever_exist",
                "no_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
