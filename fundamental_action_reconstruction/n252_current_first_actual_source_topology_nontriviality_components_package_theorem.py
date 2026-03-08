#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P232 = ROOT / "generated" / "p232_current_actual_source_topology_nontriviality_components_package_probe_summary.json"
OUT = ROOT / "generated" / "n252_current_first_actual_source_topology_nontriviality_components_package_theorem_summary.json"


def main() -> None:
    p232 = json.loads(IN_P232.read_text(encoding="utf-8"))

    passed = (
        p232["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_COMPONENTS_PACKAGE_BELOW_ACTUAL_ASSEMBLY_AND_FULL_NONTRIVIALITY_AFTER_P232"
        and p232["actual_nontriviality_components_package_exported"] is True
        and p232["actual_nontriviality_assembly_witness_exported"] is False
        and p232["actual_nonzero_flow_discharged"] is True
        and p232["barrier_protected_sign_discharged"] is True
        and p232["actual_observer_free_scope_discharged"] is True
        and p232["full_source_topology_nontriviality_discharged"] is False
        and p232["basis_independence_discharged"] is False
        and p232["qw2191_quotient_safe_discharged"] is False
        and p232["current_selector_closure"] is False
        and p232["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N252",
        "status": (
            "N252_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_COMPONENTS_PACKAGE_THEOREM_NO_FALSE_PASS"
            if passed
            else "N252_FAIL"
        ),
        "package_name": p232["package_name"],
        "refines_future_package": p232["refines_future_package"],
        "components": p232["components"],
        "observer_role": p232["observer_role"],
        "actual_nonzero_flow_discharged": p232["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": p232["barrier_protected_sign_discharged"],
        "actual_observer_free_scope_discharged": p232["actual_observer_free_scope_discharged"],
        "actual_nontriviality_components_package_exported": p232["actual_nontriviality_components_package_exported"],
        "actual_nontriviality_assembly_witness_exported": p232["actual_nontriviality_assembly_witness_exported"],
        "full_source_topology_nontriviality_discharged": p232["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p232["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p232["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p232["current_selector_closure"],
        "current_global_qw2191_discharge": p232["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
