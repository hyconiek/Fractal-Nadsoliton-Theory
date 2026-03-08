#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P218 = ROOT / "generated" / "p218_current_source_topology_nonzero_flow_subtarget_probe_summary.json"
OUT = ROOT / "generated" / "n238_current_first_source_topology_nonzero_flow_subtarget_theorem_summary.json"


def main() -> None:
    p218 = json.loads(IN_P218.read_text())

    passed = (
        p218["status"]
        == "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_NONZERO_FLOW_SUBTARGET_BELOW_ACTUAL_NONZERO_FLOW_DISCHARGE_AFTER_P218"
        and p218["future_only"] is True
        and p218["actual_nonzero_flow_discharged"] is False
        and p218["barrier_protected_sign_discharged"] is False
        and p218["full_source_topology_nontriviality_discharged"] is False
        and p218["basis_independence_discharged"] is False
        and p218["qw2191_quotient_safe_discharged"] is False
        and p218["current_selector_closure"] is False
        and p218["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N238",
        "status": (
            "N238_DISCHARGED_CURRENT_FIRST_SOURCE_TOPOLOGY_NONZERO_FLOW_SUBTARGET_THEOREM_NO_FALSE_PASS"
            if passed
            else "N238_FAIL"
        ),
        "input_packet": p218["input_packet"],
        "subtarget": p218["subtarget"],
        "codomain": p218["codomain"],
        "observer_role": p218["observer_role"],
        "future_only": p218["future_only"],
        "actual_nonzero_flow_discharged": p218["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": p218["barrier_protected_sign_discharged"],
        "full_source_topology_nontriviality_discharged": p218["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p218["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p218["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p218["current_selector_closure"],
        "current_global_qw2191_discharge": p218["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
