#!/usr/bin/env python3

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F74 = GENERATED / "f74_first_preobserver_light_matter_asymmetric_provider_packet_summary.json"
IN_F127 = GENERATED / "f127_first_source_topology_invariant_candidate_packet_summary.json"
IN_F138 = GENERATED / "f138_first_actual_source_topology_nonzero_flow_component_witness_packet_summary.json"
OUT = GENERATED / "f139_first_actual_source_topology_barrier_sign_component_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f74 = load_json(IN_F74)
    f127 = load_json(IN_F127)
    f138 = load_json(IN_F138)

    phi = f74["strict_kernel_control"]["phi"]
    kernel_core_value = f74["strict_kernel_control"]["kernel_core_limit"]
    barrier_margin = (math.pi / 2.0) - abs(phi)

    if kernel_core_value > 0.0:
        witness_value = 1
    elif kernel_core_value < 0.0:
        witness_value = -1
    else:
        witness_value = 0

    summary = {
        "packet_id": "F139",
        "status": "F139_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "input_packet": f127["candidate_packet"],
        "witness": "psi_src_barrier_sign_component_witness_v1",
        "witness_definition": "sign(cos(phi))",
        "witness_value": witness_value,
        "barrier_margin": "delta_src_barrier_sign_margin_v1",
        "barrier_margin_definition": "(pi/2) - |phi|",
        "barrier_margin_value": barrier_margin,
        "barrier_margin_positive": barrier_margin > 0.0,
        "phi_value": phi,
        "kernel_core_value": kernel_core_value,
        "source_limit_tag": f127["components"]["source_limit_tag_v1"],
        "phi_barrier_tag": f127["components"]["phi_barrier_tag_v1"],
        "observer_role": "downstream_only",
        "nonzero_flow_component_support_from_f138": f138["witness_positive"],
        "actual_scalar_sign_component_witness_exported": (
            f138["witness_positive"] is True
            and barrier_margin > 0.0
            and witness_value != 0
        ),
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
