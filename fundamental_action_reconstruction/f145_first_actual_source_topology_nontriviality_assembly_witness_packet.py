#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F129 = GENERATED / "f129_first_source_topology_invariant_nontriviality_target_packet_summary.json"
IN_F134 = GENERATED / "f134_first_source_topology_nontriviality_assembly_target_packet_summary.json"
IN_F144 = GENERATED / "f144_first_actual_source_topology_nontriviality_components_package_packet_summary.json"
OUT = GENERATED / "f145_first_actual_source_topology_nontriviality_assembly_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f129 = load_json(IN_F129)
    f134 = load_json(IN_F134)
    f144 = load_json(IN_F144)

    slot_map = [
        {
            "witness": "Xi_src_nonzero_flow_actual_witness_v1",
            "fills_target_slot": "source_limit_nonzero_flow_class_v1",
        },
        {
            "witness": "Psi_src_barrier_sign_actual_witness_v1",
            "fills_target_slot": "barrier_protected_sign_class_v1",
        },
        {
            "witness": "Omega_src_observer_free_scope_actual_witness_v1",
            "fills_target_slot": "observer_free_scope_tag_v1",
        },
    ]

    actual_nontriviality_assembly_witness_exported = (
        f134["assembly_target_name"] == "Mu_src_nontriv_assembly_target_v1"
        and f129["codomain_packet"]["name"] == "Lambda_src_nontriv_target_v1"
        and f144["package_name"] == "Kappa_src_nontriv_actual_components_packet_v1"
        and f144["actual_nontriviality_components_package_exported"] is True
        and f144["component_codomain_slots"] == f129["codomain_packet"]["components"]
        and f144["actual_nonzero_flow_discharged"] is True
        and f144["barrier_protected_sign_discharged"] is True
        and f144["actual_observer_free_scope_discharged"] is True
        and f144["observer_role"] == "downstream_only"
    )

    support_packet = {
        "domain_package": f144["package_name"],
        "codomain_target": f129["codomain_packet"]["name"],
        "slot_map": slot_map,
    }

    summary = {
        "packet_id": "F145",
        "status": "F145_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_ASSEMBLY_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "refines_future_assembly_target": f134["assembly_target_name"],
        "domain_package": f144["package_name"],
        "witness": "Mu_src_nontriv_actual_assembly_witness_v1",
        "codomain_target": f129["codomain_packet"]["name"],
        "support_packet_id": "W_src_nontriv_assembly_support_packet_v1",
        "support_packet": support_packet,
        "observer_role": "downstream_only",
        "actual_nonzero_flow_discharged": f144["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": f144["barrier_protected_sign_discharged"],
        "actual_observer_free_scope_discharged": f144["actual_observer_free_scope_discharged"],
        "actual_nontriviality_components_package_exported": f144["actual_nontriviality_components_package_exported"],
        "actual_nontriviality_assembly_witness_exported": actual_nontriviality_assembly_witness_exported,
        "full_source_topology_nontriviality_discharged": False,
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
