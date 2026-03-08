#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F146 = ROOT / "generated" / "f146_first_actual_source_topology_full_nontriviality_witness_packet_summary.json"
OUT = ROOT / "generated" / "p234_current_actual_source_topology_full_nontriviality_witness_probe_summary.json"


def main() -> None:
    f146 = json.loads(IN_F146.read_text(encoding="utf-8"))

    passed = (
        f146["input_packet"] == "tau_src_candidate_v1"
        and f146["witness"] == "Theta_src_nontriv_actual_discharge_witness_v1"
        and f146["codomain_target"] == "actual_full_source_topology_nontriviality_discharge_target_v1"
        and f146["actual_full_source_topology_nontriviality_witness_exported"] is True
        and f146["full_source_topology_nontriviality_discharged"] is True
        and f146["actual_selector_witness_exported"] is False
        and f146["observer_role"] == "downstream_only"
        and f146["basis_independence_discharged"] is False
        and f146["qw2191_quotient_safe_discharged"] is False
        and f146["current_selector_closure"] is False
        and f146["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P234",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_FULL_NONTRIVIALITY_WITNESS_BELOW_SELECTOR_WITNESS_AND_BASIS_INDEPENDENCE_AFTER_P234"
            if passed
            else "P234_FAIL"
        ),
        "input_packet": f146["input_packet"],
        "witness": f146["witness"],
        "codomain_target": f146["codomain_target"],
        "support_packet_id": f146["support_packet_id"],
        "observer_role": f146["observer_role"],
        "actual_nonzero_flow_discharged": f146["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": f146["barrier_protected_sign_discharged"],
        "actual_observer_free_scope_discharged": f146["actual_observer_free_scope_discharged"],
        "actual_nontriviality_components_package_exported": f146["actual_nontriviality_components_package_exported"],
        "actual_nontriviality_assembly_witness_exported": f146["actual_nontriviality_assembly_witness_exported"],
        "actual_full_source_topology_nontriviality_witness_exported": f146["actual_full_source_topology_nontriviality_witness_exported"],
        "full_source_topology_nontriviality_discharged": f146["full_source_topology_nontriviality_discharged"],
        "actual_selector_witness_exported": f146["actual_selector_witness_exported"],
        "basis_independence_discharged": f146["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f146["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f146["current_selector_closure"],
        "current_global_qw2191_discharge": f146["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
