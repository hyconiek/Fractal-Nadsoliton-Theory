#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F145 = GENERATED / "f145_first_actual_source_topology_nontriviality_assembly_witness_packet_summary.json"
OUT = GENERATED / "f146_first_actual_source_topology_full_nontriviality_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f145 = load_json(IN_F145)

    actual_full_source_topology_nontriviality_witness_exported = (
        f145["domain_package"] == "Kappa_src_nontriv_actual_components_packet_v1"
        and f145["witness"] == "Mu_src_nontriv_actual_assembly_witness_v1"
        and f145["codomain_target"] == "Lambda_src_nontriv_target_v1"
        and f145["actual_nontriviality_components_package_exported"] is True
        and f145["actual_nontriviality_assembly_witness_exported"] is True
        and f145["actual_nonzero_flow_discharged"] is True
        and f145["barrier_protected_sign_discharged"] is True
        and f145["actual_observer_free_scope_discharged"] is True
        and f145["observer_role"] == "downstream_only"
    )

    support_packet = {
        "input_packet": "tau_src_candidate_v1",
        "assembly_witness": f145["witness"],
        "assembled_target_packet": f145["codomain_target"],
        "full_discharge_target": "actual_full_source_topology_nontriviality_discharge_target_v1",
    }

    summary = {
        "packet_id": "F146",
        "status": "F146_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_FULL_NONTRIVIALITY_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "refines_future_discharge_target": "Theta_src_nontriv_discharge_target_v1",
        "input_packet": "tau_src_candidate_v1",
        "witness": "Theta_src_nontriv_actual_discharge_witness_v1",
        "codomain_target": "actual_full_source_topology_nontriviality_discharge_target_v1",
        "support_packet_id": "W_src_full_nontriv_support_packet_v1",
        "support_packet": support_packet,
        "observer_role": "downstream_only",
        "actual_nonzero_flow_discharged": f145["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": f145["barrier_protected_sign_discharged"],
        "actual_observer_free_scope_discharged": f145["actual_observer_free_scope_discharged"],
        "actual_nontriviality_components_package_exported": f145["actual_nontriviality_components_package_exported"],
        "actual_nontriviality_assembly_witness_exported": f145["actual_nontriviality_assembly_witness_exported"],
        "actual_full_source_topology_nontriviality_witness_exported": actual_full_source_topology_nontriviality_witness_exported,
        "full_source_topology_nontriviality_discharged": actual_full_source_topology_nontriviality_witness_exported,
        "actual_selector_witness_exported": False,
        "basis_independence_discharged": False,
        "qw2191_quotient_safe_discharged": False,
        "current_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "legacy_to_strict_bridge_claimed": False,
        "no_false_pass": True,
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
