from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text())


def main() -> None:
    c32 = load("c32_cross_pair_overlap_scalar_degeneracy_audit_summary.json")
    c33 = load("c33_local_phase_export_class_audit_summary.json")
    c34 = load("c34_local_reduced_representative_class_audit_summary.json")
    c49 = load("c49_conditional_populated_instance_schema_audit_summary.json")
    c50 = load("c50_actual_phase_source_skeleton_audit_summary.json")
    c51 = load("c51_strict_to_axiom_source_bridge_spec_audit_summary.json")
    t4 = load("t4_strict_core_export_completeness_principle_theorem_spec_summary.json")

    route_family = {
        "raw_cross_pair_overlap_route": c32["result"]["raw_cross_pair_overlap_scalar_route"],
        "local_phase_formula_route": c33["result"]["explicit_exported_theta_1_theta_2_for_actual_pair_frames"],
        "local_representative_route": c34["result"]["explicit_exported_theta_1_theta_2_for_actual_pair_frames"],
        "downstream_populated_instance_route": c49["findings"]["actual_theta_supply"],
        "strict_core_source_skeleton_route": c50["findings"]["strict_core_minimal_source_skeleton"],
        "strict_to_axiom_bridge_route": c51["findings"]["strict_to_axiom_source_bridge_spec"],
    }

    out = {
        "status": "T5_EXECUTED_T4_DISCHARGE_ATTEMPT_BLOCKED_BY_NO_ROUTE_FAMILY_CLOSURE_CERTIFICATE",
        "goal": "Attempt discharge of T4 by checking whether the current audited theta-export route family is already formally exhaustive for the present strict-core selector track.",
        "supports": {
            "C32": "raw overlap route classified",
            "C33": "formula-class route classified",
            "C34": "representative route classified",
            "C49": "downstream route classified",
            "C50": "source skeleton route classified",
            "C51": "fallback bridge route classified",
            "T4": t4["target_theorem"]
        },
        "route_family": route_family,
        "attempt_result": {
            "route_wise_classification_strength": "strong",
            "formal_exhaustiveness_certificate": "not_present",
            "route_universe_declaration": "not_present",
            "t4_discharged": False
        },
        "residual_blockers": {
            "T5_B1": "no_formal_route_family_closure_certificate_or_route_universe_declaration_showing_that_the_audited_family_C32_C33_C34_C49_C50_C51_exhausts_all_current_strict_core_actual_theta_export_routes_for_the_selector_track"
        },
        "forbidden_claims": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_t4_is_already_proved",
            "no_claim_that_the_route_family_is_already_proved_exhaustive",
            "no_claim_that_qw_2191_is_discharged"
        ]
    }

    out_path = GEN / "t5_export_completeness_principle_discharge_attempt_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(out_path)


if __name__ == "__main__":
    main()
