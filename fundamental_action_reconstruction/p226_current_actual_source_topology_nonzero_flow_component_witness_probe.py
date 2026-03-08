#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F138 = ROOT / "generated" / "f138_first_actual_source_topology_nonzero_flow_component_witness_packet_summary.json"
OUT = ROOT / "generated" / "p226_current_actual_source_topology_nonzero_flow_component_witness_probe_summary.json"


def main() -> None:
    f138 = json.loads(IN_F138.read_text())

    passed = (
        f138["input_packet"] == "tau_src_candidate_v1"
        and f138["witness"] == "xi_src_nonzero_flow_component_witness_v1"
        and f138["witness_positive"] is True
        and f138["actual_scalar_component_witness_exported"] is True
        and f138["observer_role"] == "downstream_only"
        and f138["barrier_protected_sign_discharged"] is False
        and f138["full_source_topology_nontriviality_discharged"] is False
        and f138["basis_independence_discharged"] is False
        and f138["qw2191_quotient_safe_discharged"] is False
        and f138["current_selector_closure"] is False
        and f138["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P226",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_BELOW_BARRIER_SIGN_AND_FULL_NONTRIVIALITY_AFTER_P226"
            if passed
            else "P226_FAIL"
        ),
        "input_packet": f138["input_packet"],
        "witness": f138["witness"],
        "witness_definition": f138["witness_definition"],
        "witness_value": f138["witness_value"],
        "witness_positive": f138["witness_positive"],
        "observer_role": f138["observer_role"],
        "actual_scalar_component_witness_exported": f138["actual_scalar_component_witness_exported"],
        "barrier_protected_sign_discharged": f138["barrier_protected_sign_discharged"],
        "full_source_topology_nontriviality_discharged": f138["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f138["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f138["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f138["current_selector_closure"],
        "current_global_qw2191_discharge": f138["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
