#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F144 = ROOT / "generated" / "f144_first_actual_source_topology_nontriviality_components_package_packet_summary.json"
OUT = ROOT / "generated" / "p232_current_actual_source_topology_nontriviality_components_package_probe_summary.json"


def main() -> None:
    f144 = json.loads(IN_F144.read_text(encoding="utf-8"))

    passed = (
        f144["package_name"] == "Kappa_src_nontriv_actual_components_packet_v1"
        and f144["refines_future_package"] == "Kappa_src_nontriv_components_packet_v1"
        and f144["actual_nontriviality_components_package_exported"] is True
        and f144["actual_nontriviality_assembly_witness_exported"] is False
        and f144["actual_nonzero_flow_discharged"] is True
        and f144["barrier_protected_sign_discharged"] is True
        and f144["actual_observer_free_scope_discharged"] is True
        and f144["observer_role"] == "downstream_only"
        and f144["full_source_topology_nontriviality_discharged"] is False
        and f144["basis_independence_discharged"] is False
        and f144["qw2191_quotient_safe_discharged"] is False
        and f144["current_selector_closure"] is False
        and f144["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P232",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_COMPONENTS_PACKAGE_BELOW_ACTUAL_ASSEMBLY_AND_FULL_NONTRIVIALITY_AFTER_P232"
            if passed
            else "P232_FAIL"
        ),
        "package_name": f144["package_name"],
        "refines_future_package": f144["refines_future_package"],
        "components": f144["components"],
        "observer_role": f144["observer_role"],
        "actual_nonzero_flow_discharged": f144["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": f144["barrier_protected_sign_discharged"],
        "actual_observer_free_scope_discharged": f144["actual_observer_free_scope_discharged"],
        "actual_nontriviality_components_package_exported": f144["actual_nontriviality_components_package_exported"],
        "actual_nontriviality_assembly_witness_exported": f144["actual_nontriviality_assembly_witness_exported"],
        "full_source_topology_nontriviality_discharged": f144["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f144["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f144["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f144["current_selector_closure"],
        "current_global_qw2191_discharge": f144["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
