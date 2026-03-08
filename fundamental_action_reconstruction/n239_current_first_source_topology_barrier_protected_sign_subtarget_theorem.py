#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P219 = ROOT / "generated" / "p219_current_source_topology_barrier_protected_sign_subtarget_probe_summary.json"
OUT = ROOT / "generated" / "n239_current_first_source_topology_barrier_protected_sign_subtarget_theorem_summary.json"


def main() -> None:
    p219 = json.loads(IN_P219.read_text())

    passed = (
        p219["status"]
        == "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_SUBTARGET_BELOW_ACTUAL_BARRIER_PROTECTED_SIGN_DISCHARGE_AFTER_P219"
        and p219["future_only"] is True
        and p219["actual_barrier_protected_sign_discharged"] is False
        and p219["actual_nonzero_flow_discharged"] is False
        and p219["full_source_topology_nontriviality_discharged"] is False
        and p219["basis_independence_discharged"] is False
        and p219["qw2191_quotient_safe_discharged"] is False
        and p219["current_selector_closure"] is False
        and p219["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N239",
        "status": (
            "N239_DISCHARGED_CURRENT_FIRST_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_SUBTARGET_THEOREM_NO_FALSE_PASS"
            if passed
            else "N239_FAIL"
        ),
        "input_packet": p219["input_packet"],
        "subtarget": p219["subtarget"],
        "codomain": p219["codomain"],
        "observer_role": p219["observer_role"],
        "future_only": p219["future_only"],
        "actual_barrier_protected_sign_discharged": p219["actual_barrier_protected_sign_discharged"],
        "actual_nonzero_flow_discharged": p219["actual_nonzero_flow_discharged"],
        "full_source_topology_nontriviality_discharged": p219["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p219["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p219["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p219["current_selector_closure"],
        "current_global_qw2191_discharge": p219["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
