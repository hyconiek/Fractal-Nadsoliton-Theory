#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f1_canonical_informational_fractal_substrate_parameter_packet.json"
OUT_SUMMARY = GENERATED / "f1_canonical_informational_fractal_substrate_parameter_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    tex = load_text("TOE_FINAL_DOCUMENTATION.tex")
    q1703 = load_json("material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1703_claims_vs_computation_audit.json")
    q1729 = load_text("material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW1729_NADSOLITON_KERNEL_CHARACTERISTICS_MAP.md")
    q1961 = load_text("material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW1961_NONCIRCULAR_GAMMA_Q_DERIVATION_MATRIX.md")

    df = 4.0 * math.log(2.0)
    alpha_geo = float(q1703["parameters"]["alpha_geo"])
    beta_tors = float(q1703["parameters"]["beta_tors"])
    hierarchy_ratio = alpha_geo / (2.0 * beta_tors)

    checks = [
        {
            "id": "tex_contains_D_f_4_ln_2_origin_layer_parameter",
            "actual": "D_f = 4 \\ln 2" in tex,
            "expected": True,
            "meaning": "TOE_FINAL_DOCUMENTATION.tex explicitly introduces D_f = 4 ln 2 as origin-layer fractal parameter",
        },
        {
            "id": "q1703_alpha_geo_matches_4_ln_2",
            "actual": round(alpha_geo - df, 12),
            "expected": 0.0,
            "meaning": "report_qw1703 fixes alpha_geo numerically at 4 ln 2",
        },
        {
            "id": "q1703_beta_tors_present",
            "actual": beta_tors,
            "expected": 0.01,
            "meaning": "report_qw1703 preserves beta_tors as current frozen damping parameter",
        },
        {
            "id": "q1729_mentions_alpha_geo_and_beta_tors_kernel_mapping",
            "actual": ("alpha_geo" in q1729) and ("beta_tors = 0.01" in q1729),
            "expected": True,
            "meaning": "QW-1729 maps alpha_geo and beta_tors into the nadsoliton -> kernel characteristics map",
        },
        {
            "id": "q1961_mentions_D_f_as_key_input",
            "actual": "D_f = 4 ln 2 = 2.772589" in q1961,
            "expected": True,
            "meaning": "QW-1961 explicitly uses D_f = 4 ln 2 as key input",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "F1",
        "lane": "canonical-ontology-supported-informational-fractal-substrate-parameter-layer",
        "goal": "restore_the_omitted_canonical_fractal_substrate_parameter_layer_D_f_alpha_geo_beta_tors_into_the_action_first_FIN_reconstruction_lane_without_false_strict_core_promotion",
        "canonical_provenance": {
            "tex_source": "TOE_FINAL_DOCUMENTATION.tex",
            "pdf_source": "TOE_FINAL_DOCUMENTATION 4.4.pdf",
            "repo_support_sources": [
                "report_qw1703_claims_vs_computation_audit.json",
                "RAPORT_QW1729_NADSOLITON_KERNEL_CHARACTERISTICS_MAP.md",
                "RAPORT_QW1961_NONCIRCULAR_GAMMA_Q_DERIVATION_MATRIX.md",
            ],
        },
        "canonical_structural_parameters": {
            "D_f": df,
            "alpha_geo": alpha_geo,
            "beta_tors": beta_tors,
            "alpha_geo_over_2beta_tors": hierarchy_ratio,
        },
        "current_status_classification": {
            "D_f_alpha_geo_identity": "canonical_ontology_supported_structural_parameter_identity",
            "beta_tors": "calibrated_then_frozen_damping_parameter_not_micro_derived",
            "alpha_geo_over_2beta_tors": "canonical_hierarchy_bridge_ratio_not_strict_first_principles_derivation",
        },
        "integration_obligations": {
            "A1": "must_expose_explicit_parameter_slot_or_constraint_for_D_f_alpha_geo_beta_tors_if_presented_as_informational_nadsoliton_action_first_ansatz",
            "A4": "must_expose_fractal_scaling_weight_D_f_and_damping_role_beta_tors_in_the_symbolic_coarse_graining_layer",
            "A8": "must_track_alpha_geo_over_2beta_tors_only_as_canonical_ontology_supported_hierarchy_structure_not_as_strict_derivation_of_G",
            "H_Kobs_lane": "explicit_alpha_geo_or_D_f_signal_scaling_integration_still_absent_and_remains_open",
        },
        "consistency_checks": checks,
        "anti_overclaim": {
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "strict_derivation_of_D_f": False,
            "strict_derivation_of_alpha_geo": False,
            "microscopic_derivation_of_beta_tors": False,
            "strict_derivation_of_alpha_EM_from_alpha_geo_over_2beta_tors": False,
            "strict_derivation_of_G_from_alpha_geo_over_2beta_tors": False,
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "F1",
        "status": "F1_EXECUTED_CANONICAL_INFORMATIONAL_FRACTAL_SUBSTRATE_PARAMETER_PACKET_NO_FALSE_PASS",
        "lane": artifact["lane"],
        "canonical_structural_parameters": artifact["canonical_structural_parameters"],
        "integration_obligations": artifact["integration_obligations"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
