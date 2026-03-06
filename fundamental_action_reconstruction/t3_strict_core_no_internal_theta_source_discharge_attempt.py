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
    c35 = load("c35_actual_phase_source_branch_audit_summary.json")
    c49 = load("c49_conditional_populated_instance_schema_audit_summary.json")
    c50 = load("c50_actual_phase_source_skeleton_audit_summary.json")
    c51 = load("c51_strict_to_axiom_source_bridge_spec_audit_summary.json")

    evidence = {
        "raw_overlap_route": c32["result"]["raw_cross_pair_overlap_scalar_route"],
        "c33_actual_theta_export": c33["result"]["explicit_exported_theta_1_theta_2_for_actual_pair_frames"],
        "c34_actual_theta_export": c34["result"]["explicit_exported_theta_1_theta_2_for_actual_pair_frames"],
        "c35_strict_core_source": c35["result"]["strict_core_actual_phase_source"],
        "c35_axiom_branch": c35["result"]["axiom_augmented_actual_phase_source_branch"],
        "c49_actual_theta_supply": c49["findings"]["actual_theta_supply"],
        "c50_source_skeleton": c50["findings"]["strict_core_minimal_source_skeleton"],
        "c51_bridge_spec": c51["findings"]["strict_to_axiom_source_bridge_spec"],
    }

    blockers = {
        "T3_B1": "no_formal_export_completeness_bridge_turning_the_current_not_shown_absent_fallback_only_audit_chain_into_a_theorem_level_strict_core_no_internal_theta_source_result"
    }

    out = {
        "status": "T3_EXECUTED_T1_DISCHARGE_ATTEMPT_BLOCKED_BY_NO_FORMAL_EXPORT_COMPLETENESS_BRIDGE",
        "goal": "Attempt discharge of T1 using the current strict-core selector-track audit chain without false PASS.",
        "supports": {
            "C32": "raw overlap route degenerate",
            "C33": "formula class exists, actual theta export not shown",
            "C34": "representative class exists, actual theta export not shown",
            "C35": "strict-core source absent, axiom branch present",
            "C49": "populated instance schema is downstream of actual theta values",
            "C50": "strict-core minimal source skeleton absent",
            "C51": "strict-to-axiom bridge spec absent"
        },
        "evidence": evidence,
        "attempt_result": {
            "direct_overlap_source_discharge": "blocked_by_degeneracy",
            "formula_class_source_discharge": "blocked_by_formula_only",
            "downstream_schema_as_source": "blocked_by_dependency_direction",
            "strict_core_minimal_source_skeleton": evidence["c50_source_skeleton"],
            "strict_to_axiom_bridge_spec": evidence["c51_bridge_spec"],
            "axiom_lane_only": evidence["c35_axiom_branch"]
        },
        "residual_blockers": blockers,
        "conclusion": {
            "t1_discharged": False,
            "best_honest_state": "audit_level_chain_strong_but_theorem_level_non_availability_not_yet_formally_lifted",
            "next_needed_object": "strict_core_export_completeness_principle_for_selector_track"
        },
        "forbidden_claims": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_t1_is_already_proved",
            "no_claim_that_qw_2191_is_discharged",
            "no_claim_that_axiom_lane_becomes_strict_core"
        ]
    }

    out_path = GEN / "t3_strict_core_no_internal_theta_source_discharge_attempt_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(out_path)


if __name__ == "__main__":
    main()
