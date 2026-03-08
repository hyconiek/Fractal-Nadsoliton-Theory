#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F145 = ROOT / "generated" / "f145_first_actual_source_topology_nontriviality_assembly_witness_packet_summary.json"
OUT = ROOT / "generated" / "p233_current_actual_source_topology_nontriviality_assembly_witness_probe_summary.json"


def main() -> None:
    f145 = json.loads(IN_F145.read_text(encoding="utf-8"))

    passed = (
        f145["domain_package"] == "Kappa_src_nontriv_actual_components_packet_v1"
        and f145["witness"] == "Mu_src_nontriv_actual_assembly_witness_v1"
        and f145["codomain_target"] == "Lambda_src_nontriv_target_v1"
        and f145["actual_nontriviality_components_package_exported"] is True
        and f145["actual_nontriviality_assembly_witness_exported"] is True
        and f145["observer_role"] == "downstream_only"
        and f145["full_source_topology_nontriviality_discharged"] is False
        and f145["basis_independence_discharged"] is False
        and f145["qw2191_quotient_safe_discharged"] is False
        and f145["current_selector_closure"] is False
        and f145["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P233",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_ASSEMBLY_WITNESS_BELOW_FULL_NONTRIVIALITY_AFTER_P233"
            if passed
            else "P233_FAIL"
        ),
        "domain_package": f145["domain_package"],
        "witness": f145["witness"],
        "codomain_target": f145["codomain_target"],
        "support_packet_id": f145["support_packet_id"],
        "observer_role": f145["observer_role"],
        "actual_nonzero_flow_discharged": f145["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": f145["barrier_protected_sign_discharged"],
        "actual_observer_free_scope_discharged": f145["actual_observer_free_scope_discharged"],
        "actual_nontriviality_components_package_exported": f145["actual_nontriviality_components_package_exported"],
        "actual_nontriviality_assembly_witness_exported": f145["actual_nontriviality_assembly_witness_exported"],
        "full_source_topology_nontriviality_discharged": f145["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f145["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f145["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f145["current_selector_closure"],
        "current_global_qw2191_discharge": f145["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
