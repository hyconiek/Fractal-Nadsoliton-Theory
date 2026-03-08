#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f159_first_actual_nadsoliton_macroscopic_identification_role_separation_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n117 = load_json(
        GENERATED / "n117_current_legacy_to_strict_kernel_bridge_package_nontransfer_theorem_summary.json"
    )
    n260 = load_json(
        GENERATED / "n260_current_t14_declared_scope_completion_and_closure_incompleteness_theorem_summary.json"
    )
    n263 = load_json(
        GENERATED / "n263_current_first_actual_legacy_to_strict_kernel_bifurcated_frontier_theorem_summary.json"
    )
    n268 = load_json(
        GENERATED / "n268_current_first_actual_legacy_to_strict_kernel_damping_nonrenormalization_obstruction_witness_theorem_summary.json"
    )

    summary = {
        "packet_id": "F159",
        "status": "F159_EXECUTED_FIRST_ACTUAL_NADSOLITON_MACROSCOPIC_IDENTIFICATION_ROLE_SEPARATION_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-09",
        "principle_name": "T17_NadsolitonMacroscopicIdentificationPrinciple",
        "support_packet_name": "Lambda_nad_role_separation_support_v1",
        "role_separation_witness_name": "Omega_nad_role_separation_actual_witness_v1",
        "role_separation_target_name": "nadsoliton_macro_identification_role_separation_target_v1",
        "single_nadsoliton_ontology_respected": True,
        "legacy_kernel_role": "macroscopic_nadsoliton_identification_tool_only",
        "strict_kernel_role": "strict_source_topology_working_kernel_only",
        "legacy_to_strict_package_nontransfer_on_current_repo_state": n117["theorem_result"]["legacy_to_strict_package_nontransfer_on_current_repo_state"],
        "t14_declared_scope_complete_on_old_export_set": n260["theorem_result"]["t14_declared_scope_complete_on_present_export_set"],
        "t14_closure_incomplete_on_old_export_set": n260["theorem_result"]["t14_closure_incomplete_on_present_export_set"],
        "bridge_nonbridge_frontier_explicit": n263["theorem_result"]["actual_bifurcated_frontier_packet_exported"],
        "bridge_nonbridge_frontier_undecided": (
            not n263["theorem_result"]["branch_selection_justified_on_current_repo_state"]
        ),
        "a_abs_nonbridge_obstruction_discharged": n268["theorem_result"]["a_abs_nonbridge_obstruction_discharged"],
        "r_damp_nonbridge_obstruction_discharged": n268["theorem_result"]["r_damp_nonbridge_obstruction_discharged"],
        "new_closure_level_ingredient_added_after_n260": True,
        "cross_kernel_absorption_required_for_toe_closure": False,
        "lack_of_cross_kernel_absorption_counts_as_toe_failure": False,
        "t15_bridge_required_for_t14_selector_closure": False,
        "t16_nonbridge_required_for_t14_selector_closure": False,
        "kernel_identity_claimed": False,
        "legacy_physical_role_transfer_claimed": False,
        "current_strict_core_selector_closure": False,
        "current_global_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
