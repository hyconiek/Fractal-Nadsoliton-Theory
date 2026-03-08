#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P229 = ROOT / "generated" / "p229_current_actual_source_topology_barrier_protected_sign_witness_probe_summary.json"
OUT = ROOT / "generated" / "n249_current_first_actual_source_topology_barrier_protected_sign_witness_theorem_summary.json"


def main() -> None:
    p229 = json.loads(IN_P229.read_text(encoding="utf-8"))

    passed = (
        p229["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_BELOW_FULL_SOURCE_TOPOLOGY_NONTRIVIALITY_AFTER_P229"
        and p229["actual_barrier_protected_sign_witness_exported"] is True
        and p229["barrier_protected_sign_discharged"] is True
        and p229["full_source_topology_nontriviality_discharged"] is False
        and p229["basis_independence_discharged"] is False
        and p229["qw2191_quotient_safe_discharged"] is False
        and p229["current_selector_closure"] is False
        and p229["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N249",
        "status": (
            "N249_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N249_FAIL"
        ),
        "input_packet": p229["input_packet"],
        "witness": p229["witness"],
        "codomain": p229["codomain"],
        "support_packet_id": p229["support_packet_id"],
        "observer_role": p229["observer_role"],
        "barrier_protected_sign_discharged": p229["barrier_protected_sign_discharged"],
        "full_source_topology_nontriviality_discharged": p229["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p229["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p229["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p229["current_selector_closure"],
        "current_global_qw2191_discharge": p229["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
