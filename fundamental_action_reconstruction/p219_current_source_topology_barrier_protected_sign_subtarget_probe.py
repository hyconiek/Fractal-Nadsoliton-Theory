#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F131 = ROOT / "generated" / "f131_first_source_topology_barrier_protected_sign_subtarget_packet_summary.json"
OUT = ROOT / "generated" / "p219_current_source_topology_barrier_protected_sign_subtarget_probe_summary.json"


def main() -> None:
    f131 = json.loads(IN_F131.read_text())

    passed = (
        f131["input_packet"] == "tau_src_candidate_v1"
        and f131["subtarget"] == "Psi_src_barrier_sign_target_v1"
        and f131["codomain"] == "barrier_protected_sign_class_v1"
        and f131["observer_role"] == "downstream_only"
        and f131["future_only"] is True
        and f131["actual_barrier_protected_sign_discharged"] is False
        and f131["actual_nonzero_flow_discharged"] is False
        and f131["full_source_topology_nontriviality_discharged"] is False
        and f131["basis_independence_discharged"] is False
        and f131["qw2191_quotient_safe_discharged"] is False
        and f131["current_selector_closure"] is False
        and f131["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P219",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_SUBTARGET_BELOW_ACTUAL_BARRIER_PROTECTED_SIGN_DISCHARGE_AFTER_P219"
            if passed
            else "P219_FAIL"
        ),
        "input_packet": f131["input_packet"],
        "subtarget": f131["subtarget"],
        "codomain": f131["codomain"],
        "observer_role": f131["observer_role"],
        "future_only": f131["future_only"],
        "actual_barrier_protected_sign_discharged": f131["actual_barrier_protected_sign_discharged"],
        "actual_nonzero_flow_discharged": f131["actual_nonzero_flow_discharged"],
        "full_source_topology_nontriviality_discharged": f131["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f131["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f131["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f131["current_selector_closure"],
        "current_global_qw2191_discharge": f131["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
