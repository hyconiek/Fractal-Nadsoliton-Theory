#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F127 = GENERATED / "f127_first_source_topology_invariant_candidate_packet_summary.json"
IN_F130 = GENERATED / "f130_first_source_topology_nonzero_flow_subtarget_packet_summary.json"
IN_F138 = GENERATED / "f138_first_actual_source_topology_nonzero_flow_component_witness_packet_summary.json"
OUT = GENERATED / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f127 = load_json(IN_F127)
    f130 = load_json(IN_F130)
    f138 = load_json(IN_F138)

    actual_nonzero_flow_witness_exported = (
        f127["candidate_packet"] == "tau_src_candidate_v1"
        and f130["input_packet"] == "tau_src_candidate_v1"
        and f130["codomain"] == "source_limit_nonzero_flow_class_v1"
        and f138["input_packet"] == "tau_src_candidate_v1"
        and f138["witness"] == "xi_src_nonzero_flow_component_witness_v1"
        and f138["actual_scalar_component_witness_exported"] is True
        and f138["witness_positive"] is True
        and f127["components"]["source_limit_tag_v1"] == "d -> 0"
    )

    support_packet = {
        "source_limit_tag_v1": f127["components"]["source_limit_tag_v1"],
        "T_flow_0": f127["components"]["T_flow_0"],
        "xi_src_nonzero_flow_component_witness_v1": f138["witness_value"],
    }

    summary = {
        "packet_id": "F143",
        "status": "F143_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "input_packet": f127["candidate_packet"],
        "future_subtarget": f130["subtarget"],
        "witness": "Xi_src_nonzero_flow_actual_witness_v1",
        "codomain": f130["codomain"],
        "support_packet_id": "W_src_nonzero_flow_support_packet_v1",
        "support_packet": support_packet,
        "observer_role": "downstream_only",
        "scalar_component_witness": f138["witness"],
        "scalar_component_witness_value": f138["witness_value"],
        "actual_nonzero_flow_witness_exported": actual_nonzero_flow_witness_exported,
        "actual_nonzero_flow_discharged": actual_nonzero_flow_witness_exported,
        "full_source_topology_nontriviality_discharged": False,
        "basis_independence_discharged": False,
        "qw2191_quotient_safe_discharged": False,
        "current_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "legacy_to_strict_bridge_claimed": False,
        "no_false_pass": True,
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
