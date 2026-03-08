#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F127 = GENERATED / "f127_first_source_topology_invariant_candidate_packet_summary.json"
IN_F139 = GENERATED / "f139_first_actual_source_topology_barrier_sign_component_witness_packet_summary.json"
IN_F140 = GENERATED / "f140_first_actual_source_topology_local_barrier_sign_stability_witness_packet_summary.json"
OUT = GENERATED / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f127 = load_json(IN_F127)
    f139 = load_json(IN_F139)
    f140 = load_json(IN_F140)

    actual_barrier_protected_sign_witness_exported = (
        f139["witness_value"] == 1
        and f139["barrier_margin_positive"] is True
        and f140["actual_local_sign_stability_witness_exported"] is True
        and f140["interval_inside_positive_cos_domain"] is True
    )

    support_packet = {
        "phi_barrier_tag_v1": f127["components"]["phi_barrier_tag_v1"],
        "delta_src_barrier_sign_margin_v1": f139["barrier_margin_value"],
        "epsilon_src_local_barrier_radius_v1": f140["local_radius_value"],
        "psi_src_barrier_sign_component_witness_v1": f139["witness_value"],
        "chi_src_local_barrier_sign_stability_witness_v1": f140["witness_definition"],
    }

    summary = {
        "packet_id": "F141",
        "status": "F141_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "input_packet": f127["candidate_packet"],
        "witness": "Psi_src_barrier_sign_actual_witness_v1",
        "codomain": "barrier_protected_sign_class_v1",
        "support_packet_id": "W_src_barrier_sign_support_packet_v1",
        "support_packet": support_packet,
        "observer_role": "downstream_only",
        "actual_barrier_protected_sign_witness_exported": actual_barrier_protected_sign_witness_exported,
        "barrier_protected_sign_discharged": actual_barrier_protected_sign_witness_exported,
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
