#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P230 = ROOT / "generated" / "p230_current_actual_source_topology_observer_free_scope_witness_probe_summary.json"
OUT = ROOT / "generated" / "n250_current_first_actual_source_topology_observer_free_scope_witness_theorem_summary.json"


def main() -> None:
    p230 = json.loads(IN_P230.read_text(encoding="utf-8"))

    passed = (
        p230["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_WITNESS_BELOW_FULL_SOURCE_TOPOLOGY_NONTRIVIALITY_AFTER_P230"
        and p230["actual_observer_free_scope_witness_exported"] is True
        and p230["actual_observer_free_scope_discharged"] is True
        and p230["full_source_topology_nontriviality_discharged"] is False
        and p230["basis_independence_discharged"] is False
        and p230["qw2191_quotient_safe_discharged"] is False
        and p230["current_selector_closure"] is False
        and p230["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N250",
        "status": (
            "N250_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N250_FAIL"
        ),
        "input_packet": p230["input_packet"],
        "witness": p230["witness"],
        "codomain": p230["codomain"],
        "support_packet_id": p230["support_packet_id"],
        "observer_role": p230["observer_role"],
        "actual_observer_free_scope_discharged": p230["actual_observer_free_scope_discharged"],
        "barrier_protected_sign_discharged": p230["barrier_protected_sign_discharged"],
        "actual_nonzero_flow_discharged": p230["actual_nonzero_flow_discharged"],
        "full_source_topology_nontriviality_discharged": p230["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p230["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p230["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p230["current_selector_closure"],
        "current_global_qw2191_discharge": p230["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
