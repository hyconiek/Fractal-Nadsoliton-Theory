#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F142 = ROOT / "generated" / "f142_first_actual_source_topology_observer_free_scope_witness_packet_summary.json"
OUT = ROOT / "generated" / "p230_current_actual_source_topology_observer_free_scope_witness_probe_summary.json"


def main() -> None:
    f142 = json.loads(IN_F142.read_text(encoding="utf-8"))

    passed = (
        f142["input_packet"] == "tau_src_candidate_v1"
        and f142["witness"] == "Omega_src_observer_free_scope_actual_witness_v1"
        and f142["codomain"] == "observer_free_scope_tag_v1"
        and f142["actual_observer_free_scope_witness_exported"] is True
        and f142["actual_observer_free_scope_discharged"] is True
        and f142["observer_role"] == "downstream_only"
        and f142["full_source_topology_nontriviality_discharged"] is False
        and f142["basis_independence_discharged"] is False
        and f142["qw2191_quotient_safe_discharged"] is False
        and f142["current_selector_closure"] is False
        and f142["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P230",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_WITNESS_BELOW_FULL_SOURCE_TOPOLOGY_NONTRIVIALITY_AFTER_P230"
            if passed
            else "P230_FAIL"
        ),
        "input_packet": f142["input_packet"],
        "witness": f142["witness"],
        "codomain": f142["codomain"],
        "support_packet_id": f142["support_packet_id"],
        "observer_role": f142["observer_role"],
        "actual_observer_free_scope_witness_exported": f142["actual_observer_free_scope_witness_exported"],
        "actual_observer_free_scope_discharged": f142["actual_observer_free_scope_discharged"],
        "barrier_protected_sign_discharged": f142["barrier_protected_sign_discharged"],
        "actual_nonzero_flow_discharged": f142["actual_nonzero_flow_discharged"],
        "full_source_topology_nontriviality_discharged": f142["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f142["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f142["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f142["current_selector_closure"],
        "current_global_qw2191_discharge": f142["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
