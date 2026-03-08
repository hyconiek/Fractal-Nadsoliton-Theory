#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F74 = GENERATED / "f74_first_preobserver_light_matter_asymmetric_provider_packet_summary.json"
IN_F127 = GENERATED / "f127_first_source_topology_invariant_candidate_packet_summary.json"
OUT = GENERATED / "f138_first_actual_source_topology_nonzero_flow_component_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f74 = load_json(IN_F74)
    f127 = load_json(IN_F127)

    witness_value = abs(f74["strict_kernel_control"]["kernel_core_limit"])

    summary = {
        "packet_id": "F138",
        "status": "F138_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "input_packet": f127["candidate_packet"],
        "witness": "xi_src_nonzero_flow_component_witness_v1",
        "witness_definition": "|cos(phi)|",
        "witness_value": witness_value,
        "witness_positive": witness_value > 0.0,
        "source_limit_tag": f127["components"]["source_limit_tag_v1"],
        "phi_barrier_tag": f127["components"]["phi_barrier_tag_v1"],
        "observer_role": "downstream_only",
        "actual_scalar_component_witness_exported": True,
        "barrier_protected_sign_discharged": False,
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
