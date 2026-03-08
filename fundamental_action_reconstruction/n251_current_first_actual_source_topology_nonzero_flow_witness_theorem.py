#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P231 = ROOT / "generated" / "p231_current_actual_source_topology_nonzero_flow_witness_probe_summary.json"
OUT = ROOT / "generated" / "n251_current_first_actual_source_topology_nonzero_flow_witness_theorem_summary.json"


def main() -> None:
    p231 = json.loads(IN_P231.read_text(encoding="utf-8"))

    passed = (
        p231["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_WITNESS_BELOW_FULL_SOURCE_TOPOLOGY_NONTRIVIALITY_AFTER_P231"
        and p231["actual_nonzero_flow_witness_exported"] is True
        and p231["actual_nonzero_flow_discharged"] is True
        and p231["full_source_topology_nontriviality_discharged"] is False
        and p231["basis_independence_discharged"] is False
        and p231["qw2191_quotient_safe_discharged"] is False
        and p231["current_selector_closure"] is False
        and p231["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N251",
        "status": (
            "N251_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N251_FAIL"
        ),
        "input_packet": p231["input_packet"],
        "witness": p231["witness"],
        "codomain": p231["codomain"],
        "support_packet_id": p231["support_packet_id"],
        "observer_role": p231["observer_role"],
        "actual_nonzero_flow_discharged": p231["actual_nonzero_flow_discharged"],
        "full_source_topology_nontriviality_discharged": p231["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p231["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p231["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p231["current_selector_closure"],
        "current_global_qw2191_discharge": p231["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
