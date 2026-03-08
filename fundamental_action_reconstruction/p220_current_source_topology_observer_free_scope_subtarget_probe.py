#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F132 = ROOT / "generated" / "f132_first_source_topology_observer_free_scope_subtarget_packet_summary.json"
OUT = ROOT / "generated" / "p220_current_source_topology_observer_free_scope_subtarget_probe_summary.json"


def main() -> None:
    f132 = json.loads(IN_F132.read_text())

    passed = (
        f132["input_packet"] == "tau_src_candidate_v1"
        and f132["subtarget"] == "Omega_src_observer_free_scope_target_v1"
        and f132["codomain"] == "observer_free_scope_tag_v1"
        and f132["observer_role"] == "downstream_only"
        and f132["future_only"] is True
        and f132["actual_observer_free_scope_discharged"] is False
        and f132["actual_barrier_protected_sign_discharged"] is False
        and f132["actual_nonzero_flow_discharged"] is False
        and f132["full_source_topology_nontriviality_discharged"] is False
        and f132["basis_independence_discharged"] is False
        and f132["qw2191_quotient_safe_discharged"] is False
        and f132["current_selector_closure"] is False
        and f132["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P220",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_SUBTARGET_BELOW_ACTUAL_OBSERVER_FREE_SCOPE_DISCHARGE_AFTER_P220"
            if passed
            else "P220_FAIL"
        ),
        "input_packet": f132["input_packet"],
        "subtarget": f132["subtarget"],
        "codomain": f132["codomain"],
        "observer_role": f132["observer_role"],
        "future_only": f132["future_only"],
        "actual_observer_free_scope_discharged": f132["actual_observer_free_scope_discharged"],
        "actual_barrier_protected_sign_discharged": f132["actual_barrier_protected_sign_discharged"],
        "actual_nonzero_flow_discharged": f132["actual_nonzero_flow_discharged"],
        "full_source_topology_nontriviality_discharged": f132["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f132["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f132["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f132["current_selector_closure"],
        "current_global_qw2191_discharge": f132["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
