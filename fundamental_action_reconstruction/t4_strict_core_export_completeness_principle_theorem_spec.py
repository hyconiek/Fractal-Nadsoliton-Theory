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
    t3 = load("t3_strict_core_no_internal_theta_source_discharge_attempt_summary.json")

    out = {
        "status": "T4_PACKET_READY_STRICT_CORE_EXPORT_COMPLETENESS_PRINCIPLE_THEOREM_SPEC_NO_FALSE_PASS",
        "goal": "Write a packet-ready theorem spec for the missing export-completeness principle isolated by T3.",
        "supports": {
            "C32": c32["result"]["raw_cross_pair_overlap_scalar_route"],
            "C33": c33["result"]["explicit_exported_theta_1_theta_2_for_actual_pair_frames"],
            "C34": c34["result"]["explicit_exported_theta_1_theta_2_for_actual_pair_frames"],
            "C49": c49["findings"]["actual_theta_supply"],
            "C50": c50["findings"]["strict_core_minimal_source_skeleton"],
            "C51": c51["findings"]["strict_to_axiom_source_bridge_spec"],
            "T3": t3["residual_blockers"]["T3_B1"]
        },
        "route_family_under_test": [
            "raw_cross_pair_overlap_route",
            "formula_class_route",
            "representative_class_route",
            "downstream_populated_instance_schema",
            "strict_core_minimal_source_skeleton",
            "strict_to_axiom_fallback_bridge"
        ],
        "target_theorem": "strict_core_export_completeness_principle_for_current_selector_track",
        "residual_blockers": {
            "T4_B1": "the_export_completeness_principle_is_specified_but_not_discharged_for_the_current_strict_core_selector_track"
        },
        "forbidden_claims": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_export_completeness_is_already_proved",
            "no_claim_that_t1_is_already_discharged",
            "no_claim_that_qw_2191_is_discharged"
        ]
    }

    out_path = GEN / "t4_strict_core_export_completeness_principle_theorem_spec_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(out_path)


if __name__ == "__main__":
    main()
