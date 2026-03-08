#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F130 = ROOT / "generated" / "f130_first_source_topology_nonzero_flow_subtarget_packet_summary.json"
OUT = ROOT / "generated" / "p218_current_source_topology_nonzero_flow_subtarget_probe_summary.json"


def main() -> None:
    f130 = json.loads(IN_F130.read_text())

    passed = (
        f130["input_packet"] == "tau_src_candidate_v1"
        and f130["subtarget"] == "Xi_src_nonzero_flow_target_v1"
        and f130["codomain"] == "source_limit_nonzero_flow_class_v1"
        and f130["observer_role"] == "downstream_only"
        and f130["future_only"] is True
        and f130["actual_nonzero_flow_discharged"] is False
        and f130["barrier_protected_sign_discharged"] is False
        and f130["full_source_topology_nontriviality_discharged"] is False
        and f130["basis_independence_discharged"] is False
        and f130["qw2191_quotient_safe_discharged"] is False
        and f130["current_selector_closure"] is False
        and f130["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P218",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_NONZERO_FLOW_SUBTARGET_BELOW_ACTUAL_NONZERO_FLOW_DISCHARGE_AFTER_P218"
            if passed
            else "P218_FAIL"
        ),
        "input_packet": f130["input_packet"],
        "subtarget": f130["subtarget"],
        "codomain": f130["codomain"],
        "observer_role": f130["observer_role"],
        "future_only": f130["future_only"],
        "actual_nonzero_flow_discharged": f130["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": f130["barrier_protected_sign_discharged"],
        "full_source_topology_nontriviality_discharged": f130["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f130["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f130["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f130["current_selector_closure"],
        "current_global_qw2191_discharge": f130["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
