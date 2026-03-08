#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P220 = ROOT / "generated" / "p220_current_source_topology_observer_free_scope_subtarget_probe_summary.json"
OUT = ROOT / "generated" / "n240_current_first_source_topology_observer_free_scope_subtarget_theorem_summary.json"


def main() -> None:
    p220 = json.loads(IN_P220.read_text())

    passed = (
        p220["status"]
        == "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_SUBTARGET_BELOW_ACTUAL_OBSERVER_FREE_SCOPE_DISCHARGE_AFTER_P220"
        and p220["future_only"] is True
        and p220["actual_observer_free_scope_discharged"] is False
        and p220["actual_barrier_protected_sign_discharged"] is False
        and p220["actual_nonzero_flow_discharged"] is False
        and p220["full_source_topology_nontriviality_discharged"] is False
        and p220["basis_independence_discharged"] is False
        and p220["qw2191_quotient_safe_discharged"] is False
        and p220["current_selector_closure"] is False
        and p220["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N240",
        "status": (
            "N240_DISCHARGED_CURRENT_FIRST_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_SUBTARGET_THEOREM_NO_FALSE_PASS"
            if passed
            else "N240_FAIL"
        ),
        "input_packet": p220["input_packet"],
        "subtarget": p220["subtarget"],
        "codomain": p220["codomain"],
        "observer_role": p220["observer_role"],
        "future_only": p220["future_only"],
        "actual_observer_free_scope_discharged": p220["actual_observer_free_scope_discharged"],
        "actual_barrier_protected_sign_discharged": p220["actual_barrier_protected_sign_discharged"],
        "actual_nonzero_flow_discharged": p220["actual_nonzero_flow_discharged"],
        "full_source_topology_nontriviality_discharged": p220["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p220["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p220["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p220["current_selector_closure"],
        "current_global_qw2191_discharge": p220["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
