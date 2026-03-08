#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P226 = ROOT / "generated" / "p226_current_actual_source_topology_nonzero_flow_component_witness_probe_summary.json"
OUT = ROOT / "generated" / "n246_current_first_actual_source_topology_nonzero_flow_component_witness_theorem_summary.json"


def main() -> None:
    p226 = json.loads(IN_P226.read_text())

    passed = (
        p226["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_BELOW_BARRIER_SIGN_AND_FULL_NONTRIVIALITY_AFTER_P226"
        and p226["witness_positive"] is True
        and p226["barrier_protected_sign_discharged"] is False
        and p226["full_source_topology_nontriviality_discharged"] is False
        and p226["basis_independence_discharged"] is False
        and p226["qw2191_quotient_safe_discharged"] is False
        and p226["current_selector_closure"] is False
        and p226["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N246",
        "status": (
            "N246_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N246_FAIL"
        ),
        "input_packet": p226["input_packet"],
        "witness": p226["witness"],
        "witness_definition": p226["witness_definition"],
        "witness_value": p226["witness_value"],
        "witness_positive": p226["witness_positive"],
        "observer_role": p226["observer_role"],
        "barrier_protected_sign_discharged": p226["barrier_protected_sign_discharged"],
        "full_source_topology_nontriviality_discharged": p226["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p226["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p226["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p226["current_selector_closure"],
        "current_global_qw2191_discharge": p226["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
