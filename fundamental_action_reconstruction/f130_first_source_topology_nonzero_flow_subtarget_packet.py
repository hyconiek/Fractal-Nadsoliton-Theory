#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "f130_first_source_topology_nonzero_flow_subtarget_packet_summary.json"


def main() -> None:
    summary = {
        "packet_id": "F130",
        "status": "F130_EXECUTED_FIRST_SOURCE_TOPOLOGY_NONZERO_FLOW_SUBTARGET_PACKET_FUTURE_ROUTE_ONLY_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "input_packet": "tau_src_candidate_v1",
        "subtarget": "Xi_src_nonzero_flow_target_v1",
        "codomain": "source_limit_nonzero_flow_class_v1",
        "observer_role": "downstream_only",
        "future_only": True,
        "actual_nonzero_flow_discharged": False,
        "barrier_protected_sign_discharged": False,
        "full_source_topology_nontriviality_discharged": False,
        "basis_independence_discharged": False,
        "qw2191_quotient_safe_discharged": False,
        "current_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "legacy_to_strict_bridge_claimed": False,
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
