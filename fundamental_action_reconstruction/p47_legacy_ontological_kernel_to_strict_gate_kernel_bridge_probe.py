#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p47_legacy_ontological_kernel_to_strict_gate_kernel_bridge_probe.json"
OUT_SUMMARY = (
    GENERATED
    / "p47_legacy_ontological_kernel_to_strict_gate_kernel_bridge_probe_summary.json"
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
    legacy_tex = load_text("TOE_FINAL_DOCUMENTATION_RELEASE_4_4_LEGACY_FULL.tex")
    strict_textbook = load_text("RELEASE_4_9_TEXTBOOK_EN_PL.md")
    qw2005 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2005_PRE1700_TEX_REVALIDATION_MATRIX_EN_PL.md"
    )
    qw2049 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.md"
    )

    legacy_kernel_present = contains_all(
        legacy_tex,
        [
            "K(d) = \\alphageo",
            "\\cos(\\omega d + \\phi)",
            "1 + \\betators d",
        ],
    )
    strict_kernel_present = contains_all(
        strict_textbook,
        [
            "K(d)=\\frac{\\cos(\\omega d+\\phi)}{1+\\beta d^{\\eta}}",
            "0.18575",
            "0.16250",
            "1.80000",
        ],
    )
    strict_gate_selected_tuple_present = contains_all(
        qw2049,
        ["0.185750 / 0.162500 / 1.000000 / 1.800000"]
    )
    canonical_legacy_parameter_layer_restored = contains_all(
        json.dumps(f1, sort_keys=True),
        ["2.772588722239781", "0.01", "138.62943611198907"],
    )
    legacy_revalidation_downgrades_present = contains_all(
        qw2005,
        [
            "PARTIAL / CONFLICTED",
            "HEURISTIC / NOT STRICTLY DERIVED",
            "PARTIAL / MODEL-LEVEL",
            "MODEL-CONSISTENT, NOT FULL INDEPENDENT PROOF",
        ],
    )
    partial_retained_vs_replaced_partition_present = contains_all(
        qw2005,
        [
            "Retained:",
            "Replaced / upgraded:",
            "oscillatory+damped kernel structure idea",
            "strict closure stack",
        ],
    )

    authoritative_shared_parameter_bridge_present = (
        contains_all(strict_textbook, ["alpha_geo", "0.18575", "1.80000"])
        or contains_all(qw2049, ["alpha_geo", "1.800000"])
        or contains_all(legacy_tex, ["0.18575", "1.80000"])
    )
    explicit_beta_tors_to_beta_eta_translation_present = (
        contains_all(strict_textbook, ["beta_tors", "beta", "eta"])
        or contains_all(qw2049, ["beta_tors", "beta", "eta"])
        or contains_all(legacy_tex, ["beta_tors", "beta", "eta"])
    )
    explicit_phase_frequency_bridge_present = (
        contains_all(strict_textbook, ["pi/4", "0.18575"])
        or contains_all(strict_textbook, ["pi/6", "0.16250"])
        or contains_all(qw2049, ["pi/4", "0.185750"])
        or contains_all(qw2049, ["pi/6", "0.162500"])
        or contains_all(legacy_tex, ["0.18575", "1.80000"])
    )
    explicit_kernel_level_ontological_anchoring_witness_present = (
        contains_all(strict_textbook, ["the nadsoliton itself", "K(d)"])
        or contains_all(qw2049, ["nadsoliton", "0.185750", "1.800000"])
    )

    remaining_missing = [
        "explicit_amplitude_normalization_or_absorption_map_explaining_the_loss_of_visible_alpha_geo_between_K_legacy_ont_and_K_strict_gate",
        "explicit_renormalized_damping_map_beta_tors_to_beta_eta_or_equivalent_strict_translation_rule",
        "explicit_phase_frequency_bridge_from_pi_over_4_pi_over_6_to_0p18575_0p16250",
        "explicit_kernel_specific_retained_vs_replaced_partition_theorem_for_legacy_vs_strict_kernel_roles",
        "explicit_kernel_level_ontological_anchoring_witness_showing_that_K_strict_gate_is_an_internal_route_of_one_informational_nadsoliton",
    ]

    checks_spec = [
        {
            "id": "legacy_kernel_present_in_authoritative_legacy_source",
            "actual": legacy_kernel_present,
            "expected": True,
            "meaning": "the legacy kernel is explicitly present in canonical legacy sources",
        },
        {
            "id": "strict_kernel_present_in_release_4_9_source",
            "actual": strict_kernel_present,
            "expected": True,
            "meaning": "the strict gate working kernel is explicitly present in the Release 4.9 textbook source",
        },
        {
            "id": "strict_qw2049_selected_tuple_present",
            "actual": strict_gate_selected_tuple_present,
            "expected": True,
            "meaning": "QW-2049 exports the selected strict gate tuple",
        },
        {
            "id": "f1_restores_legacy_parameter_layer",
            "actual": canonical_legacy_parameter_layer_restored,
            "expected": True,
            "meaning": "F1 restores the canonical legacy parameter layer D_f/alpha_geo/beta_tors in FAR",
        },
        {
            "id": "q2005_downgrades_old_kernel_edifice_from_strict_proof",
            "actual": legacy_revalidation_downgrades_present,
            "expected": True,
            "meaning": "QW-2005 already downgrades several legacy formulas from strong strict-proof status",
        },
        {
            "id": "q2005_contains_only_partial_retained_vs_replaced_partition",
            "actual": partial_retained_vs_replaced_partition_present,
            "expected": True,
            "meaning": "QW-2005 contains a broad retained-vs-replaced partition, but not a kernel-specific bridge theorem",
        },
        {
            "id": "authoritative_shared_parameter_bridge_present",
            "actual": authoritative_shared_parameter_bridge_present,
            "expected": False,
            "meaning": "no authoritative source currently co-presents the legacy and strict parameter sets as one explicit bridge object",
        },
        {
            "id": "explicit_beta_tors_to_beta_eta_translation_present",
            "actual": explicit_beta_tors_to_beta_eta_translation_present,
            "expected": False,
            "meaning": "no authoritative source currently exports an explicit translation from beta_tors to beta and eta",
        },
        {
            "id": "explicit_phase_frequency_bridge_present",
            "actual": explicit_phase_frequency_bridge_present,
            "expected": False,
            "meaning": "no authoritative source currently exports an explicit bridge from the legacy phase/frequency tuple to the strict tuple",
        },
        {
            "id": "explicit_kernel_level_ontological_anchoring_witness_present",
            "actual": explicit_kernel_level_ontological_anchoring_witness_present,
            "expected": False,
            "meaning": "the strict-side authoritative sources still do not explicitly re-anchor the strict kernel as an internal route of one informational nadsoliton",
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
        "stage": "P47",
        "lane": "legacy_ontological_kernel_to_strict_gate_kernel_bridge_probe_current_repo_state_only",
        "goal": "test_whether_the_current_repo_already_exports_a_rigorous_bridge_identifying_the_legacy_ontological_kernel_with_the_later_strict_gate_kernel",
        "status": "CURRENT_REPO_EXPORTS_LEGACY_AND_STRICT_KERNELS_BUT_NO_RIGOROUS_LEGACY_TO_STRICT_KERNEL_BRIDGE_AFTER_P47",
        "reason": "the current repo explicitly exports the legacy ontological/effective kernel in canonical legacy sources, the strict gate working kernel in Release 4.9/QW-2049 sources, the restored legacy parameter layer D_f/alpha_geo/beta_tors in F1, and a partial old-claim downgrade in QW-2005, but it still does not export an explicit amplitude absorption map for alpha_geo, an explicit beta_tors to beta/eta translation rule, an explicit phase/frequency bridge, an explicit kernel-specific retained-vs-replaced theorem, or an explicit kernel-level ontological anchoring witness for the strict gate kernel",
        "kernel_objects_present": {
            "legacy_ontological_effective_kernel_present": legacy_kernel_present,
            "strict_gate_working_kernel_present": strict_kernel_present,
            "strict_gate_selected_tuple_present": strict_gate_selected_tuple_present,
            "canonical_legacy_parameter_layer_restored_in_far": canonical_legacy_parameter_layer_restored,
            "partial_legacy_revalidation_downgrade_present": legacy_revalidation_downgrades_present,
            "partial_retained_vs_replaced_partition_present": partial_retained_vs_replaced_partition_present,
        },
        "bridge_state": {
            "authoritative_shared_parameter_bridge_present": authoritative_shared_parameter_bridge_present,
            "explicit_beta_tors_to_beta_eta_translation_present": explicit_beta_tors_to_beta_eta_translation_present,
            "explicit_phase_frequency_bridge_present": explicit_phase_frequency_bridge_present,
            "explicit_kernel_level_ontological_anchoring_witness_present": explicit_kernel_level_ontological_anchoring_witness_present,
        },
        "remaining_missing_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P47",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "kernel_objects_present": artifact["kernel_objects_present"],
        "bridge_state": artifact["bridge_state"],
        "remaining_missing_objects": remaining_missing,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
