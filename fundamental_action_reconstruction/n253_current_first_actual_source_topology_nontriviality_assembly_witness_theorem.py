#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P233 = ROOT / "generated" / "p233_current_actual_source_topology_nontriviality_assembly_witness_probe_summary.json"
OUT = ROOT / "generated" / "n253_current_first_actual_source_topology_nontriviality_assembly_witness_theorem_summary.json"


def main() -> None:
    p233 = json.loads(IN_P233.read_text(encoding="utf-8"))

    passed = (
        p233["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_ASSEMBLY_WITNESS_BELOW_FULL_NONTRIVIALITY_AFTER_P233"
        and p233["actual_nontriviality_components_package_exported"] is True
        and p233["actual_nontriviality_assembly_witness_exported"] is True
        and p233["actual_nonzero_flow_discharged"] is True
        and p233["barrier_protected_sign_discharged"] is True
        and p233["actual_observer_free_scope_discharged"] is True
        and p233["full_source_topology_nontriviality_discharged"] is False
        and p233["basis_independence_discharged"] is False
        and p233["qw2191_quotient_safe_discharged"] is False
        and p233["current_selector_closure"] is False
        and p233["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N253",
        "status": (
            "N253_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_ASSEMBLY_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N253_FAIL"
        ),
        "domain_package": p233["domain_package"],
        "witness": p233["witness"],
        "codomain_target": p233["codomain_target"],
        "support_packet_id": p233["support_packet_id"],
        "observer_role": p233["observer_role"],
        "actual_nonzero_flow_discharged": p233["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": p233["barrier_protected_sign_discharged"],
        "actual_observer_free_scope_discharged": p233["actual_observer_free_scope_discharged"],
        "actual_nontriviality_components_package_exported": p233["actual_nontriviality_components_package_exported"],
        "actual_nontriviality_assembly_witness_exported": p233["actual_nontriviality_assembly_witness_exported"],
        "full_source_topology_nontriviality_discharged": p233["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p233["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p233["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p233["current_selector_closure"],
        "current_global_qw2191_discharge": p233["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
