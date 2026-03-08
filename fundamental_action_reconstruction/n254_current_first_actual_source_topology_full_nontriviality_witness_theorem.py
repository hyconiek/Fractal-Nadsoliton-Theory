#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P234 = ROOT / "generated" / "p234_current_actual_source_topology_full_nontriviality_witness_probe_summary.json"
OUT = ROOT / "generated" / "n254_current_first_actual_source_topology_full_nontriviality_witness_theorem_summary.json"


def main() -> None:
    p234 = json.loads(IN_P234.read_text(encoding="utf-8"))

    passed = (
        p234["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_FULL_NONTRIVIALITY_WITNESS_BELOW_SELECTOR_WITNESS_AND_BASIS_INDEPENDENCE_AFTER_P234"
        and p234["actual_full_source_topology_nontriviality_witness_exported"] is True
        and p234["full_source_topology_nontriviality_discharged"] is True
        and p234["actual_selector_witness_exported"] is False
        and p234["basis_independence_discharged"] is False
        and p234["qw2191_quotient_safe_discharged"] is False
        and p234["current_selector_closure"] is False
        and p234["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N254",
        "status": (
            "N254_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_FULL_NONTRIVIALITY_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N254_FAIL"
        ),
        "input_packet": p234["input_packet"],
        "witness": p234["witness"],
        "codomain_target": p234["codomain_target"],
        "support_packet_id": p234["support_packet_id"],
        "observer_role": p234["observer_role"],
        "actual_nonzero_flow_discharged": p234["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": p234["barrier_protected_sign_discharged"],
        "actual_observer_free_scope_discharged": p234["actual_observer_free_scope_discharged"],
        "actual_nontriviality_components_package_exported": p234["actual_nontriviality_components_package_exported"],
        "actual_nontriviality_assembly_witness_exported": p234["actual_nontriviality_assembly_witness_exported"],
        "actual_full_source_topology_nontriviality_witness_exported": p234["actual_full_source_topology_nontriviality_witness_exported"],
        "full_source_topology_nontriviality_discharged": p234["full_source_topology_nontriviality_discharged"],
        "actual_selector_witness_exported": p234["actual_selector_witness_exported"],
        "basis_independence_discharged": p234["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p234["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p234["current_selector_closure"],
        "current_global_qw2191_discharge": p234["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
