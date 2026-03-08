#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F143 = ROOT / "generated" / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"
OUT = ROOT / "generated" / "p231_current_actual_source_topology_nonzero_flow_witness_probe_summary.json"


def main() -> None:
    f143 = json.loads(IN_F143.read_text(encoding="utf-8"))

    passed = (
        f143["input_packet"] == "tau_src_candidate_v1"
        and f143["witness"] == "Xi_src_nonzero_flow_actual_witness_v1"
        and f143["codomain"] == "source_limit_nonzero_flow_class_v1"
        and f143["actual_nonzero_flow_witness_exported"] is True
        and f143["actual_nonzero_flow_discharged"] is True
        and f143["observer_role"] == "downstream_only"
        and f143["full_source_topology_nontriviality_discharged"] is False
        and f143["basis_independence_discharged"] is False
        and f143["qw2191_quotient_safe_discharged"] is False
        and f143["current_selector_closure"] is False
        and f143["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P231",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_WITNESS_BELOW_FULL_SOURCE_TOPOLOGY_NONTRIVIALITY_AFTER_P231"
            if passed
            else "P231_FAIL"
        ),
        "input_packet": f143["input_packet"],
        "witness": f143["witness"],
        "codomain": f143["codomain"],
        "support_packet_id": f143["support_packet_id"],
        "observer_role": f143["observer_role"],
        "actual_nonzero_flow_witness_exported": f143["actual_nonzero_flow_witness_exported"],
        "actual_nonzero_flow_discharged": f143["actual_nonzero_flow_discharged"],
        "full_source_topology_nontriviality_discharged": f143["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f143["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f143["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f143["current_selector_closure"],
        "current_global_qw2191_discharge": f143["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
