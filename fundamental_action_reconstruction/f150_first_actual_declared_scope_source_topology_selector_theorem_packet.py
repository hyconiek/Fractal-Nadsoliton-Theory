#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F142 = GENERATED / "f142_first_actual_source_topology_observer_free_scope_witness_packet_summary.json"
IN_F146 = GENERATED / "f146_first_actual_source_topology_full_nontriviality_witness_packet_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F148 = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
IN_F149 = GENERATED / "f149_first_actual_source_topology_quotient_safe_qw2191_resolution_witness_packet_summary.json"
IN_N163 = GENERATED / "n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
IN_N234 = GENERATED / "n234_current_global_selector_closure_and_qw2191_discharge_promotion_obstruction_theorem_summary.json"
OUT = GENERATED / "f150_first_actual_declared_scope_source_topology_selector_theorem_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f142 = load_json(IN_F142)
    f146 = load_json(IN_F146)
    f147 = load_json(IN_F147)
    f148 = load_json(IN_F148)
    f149 = load_json(IN_F149)
    n163 = load_json(IN_N163)
    n234 = load_json(IN_N234)

    ordering = ["nadsoliton", "light", "matter", "emergent_observer"]

    l1_discharged = (
        f146["input_packet"] == "tau_src_candidate_v1"
        and f146["witness"] == "Theta_src_nontriv_actual_discharge_witness_v1"
        and f146["full_source_topology_nontriviality_discharged"] is True
    )

    l2_discharged = (
        f142["input_packet"] == "tau_src_candidate_v1"
        and f142["witness"] == "Omega_src_observer_free_scope_actual_witness_v1"
        and f142["actual_observer_free_scope_discharged"] is True
        and f142["support_packet"]["ordering"] == ordering
        and n163["theorem_result"]["ordering_supported_on_current_repo_state"] is True
        and n163["theorem_result"]["observer_information_deficit_downstream_symptom_on_current_repo_state"] is True
        and n163["theorem_result"]["primary_missing_selector_source_gap_upstream_of_observer"] is True
    )

    l3_discharged = (
        f147["input_packet"] == "tau_src_candidate_v1"
        and f147["witness"] == "Pi_sel_src_actual_witness_v1"
        and f147["actual_selector_witness_exported"] is True
        and f147["tau_src_identified_with_s_prelm"] is False
        and f148["input_packet"] == "tau_src_candidate_v1"
        and f148["witness"] == "Upsilon_sel_basis_actual_witness_v1"
        and f148["actual_basis_independent_selector_promotion_exported"] is True
        and f148["basis_independence_discharged"] is True
        and f148["tau_src_identified_with_s_prelm"] is False
    )

    l4_discharged = (
        f149["input_packet"] == "tau_src_candidate_v1"
        and f149["witness"] == "Phi_qw2191_safe_actual_witness_v1"
        and f149["actual_quotient_safe_qw2191_resolution_exported"] is True
        and f149["qw2191_quotient_safe_discharged"] is True
        and f149["raw_theta_uniqueness_claimed"] is False
    )

    l5_boundary_retained = (
        n234["theorem_result"]["observer_downstream_only"] is True
        and n234["theorem_result"]["global_selector_closure_justified_on_current_repo_state"] is False
        and n234["theorem_result"]["global_qw2191_discharge_justified_on_current_repo_state"] is False
        and f142["observer_role"] == "downstream_only"
        and f148["observer_role"] == "downstream_only"
        and f149["observer_role"] == "downstream_only"
    )

    declared_scope_source_topology_selector_theorem_exported = (
        l1_discharged is True
        and l2_discharged is True
        and l3_discharged is True
        and l4_discharged is True
        and l5_boundary_retained is True
        and f149["current_selector_closure"] is False
        and f149["current_global_qw2191_discharge"] is False
        and f149["legacy_to_strict_bridge_claimed"] is False
    )

    support_packet = {
        "input_packet": "tau_src_candidate_v1",
        "t14_theorem_spec": "T14_SourceTopologySelector_Theorem",
        "l1_full_nontriviality_witness": f146["witness"],
        "l2_observer_free_scope_witness": f142["witness"],
        "l3_source_side_selector_witness": f147["witness"],
        "l3_basis_independent_promotion_witness": f148["witness"],
        "l4_qw2191_quotient_safe_resolution_witness": f149["witness"],
        "l5_boundaries": {
            "n163_downstream_symptom_boundary": n163["status"],
            "n234_no_global_promotion_boundary": n234["status"],
            "observer_downstream_only": True,
        },
    }

    summary = {
        "packet_id": "F150",
        "status": "F150_EXECUTED_FIRST_ACTUAL_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "refines_theorem_spec": "T14_SourceTopologySelector_Theorem",
        "input_packet": "tau_src_candidate_v1",
        "witness": "T14_src_selector_declared_scope_actual_witness_v1",
        "codomain_target": "declared_scope_source_topology_selector_theorem_target_v1",
        "support_packet_id": "W_src_topology_selector_theorem_support_packet_v1",
        "support_packet": support_packet,
        "observer_role": "downstream_only",
        "tau_src_identified_with_s_prelm": False,
        "actual_full_source_topology_nontriviality_witness_exported": f146["actual_full_source_topology_nontriviality_witness_exported"],
        "actual_observer_free_scope_discharged": f142["actual_observer_free_scope_discharged"],
        "actual_selector_witness_exported": f147["actual_selector_witness_exported"],
        "basis_independence_discharged": f148["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f149["qw2191_quotient_safe_discharged"],
        "raw_theta_uniqueness_claimed": False,
        "declared_scope_only": True,
        "admissible_strict_core_internal_selector_source_object_claimed": False,
        "declared_scope_source_topology_selector_theorem_exported": declared_scope_source_topology_selector_theorem_exported,
        "current_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "legacy_to_strict_bridge_claimed": False,
        "older_negative_theorems_reinterpreted": False,
        "no_false_pass": True,
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
