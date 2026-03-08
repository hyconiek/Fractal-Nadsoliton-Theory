#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F141 = GENERATED / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"
IN_F142 = GENERATED / "f142_first_actual_source_topology_observer_free_scope_witness_packet_summary.json"
IN_F143 = GENERATED / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"
OUT = GENERATED / "f144_first_actual_source_topology_nontriviality_components_package_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f141 = load_json(IN_F141)
    f142 = load_json(IN_F142)
    f143 = load_json(IN_F143)

    actual_nontriviality_components_package_exported = (
        f143["actual_nonzero_flow_discharged"] is True
        and f141["barrier_protected_sign_discharged"] is True
        and f142["actual_observer_free_scope_discharged"] is True
        and f143["observer_role"] == "downstream_only"
        and f141["observer_role"] == "downstream_only"
        and f142["observer_role"] == "downstream_only"
    )

    summary = {
        "packet_id": "F144",
        "status": "F144_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_COMPONENTS_PACKAGE_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "package_name": "Kappa_src_nontriv_actual_components_packet_v1",
        "refines_future_package": "Kappa_src_nontriv_components_packet_v1",
        "components": [
            "Xi_src_nonzero_flow_actual_witness_v1",
            "Psi_src_barrier_sign_actual_witness_v1",
            "Omega_src_observer_free_scope_actual_witness_v1",
        ],
        "component_codomain_slots": [
            "source_limit_nonzero_flow_class_v1",
            "barrier_protected_sign_class_v1",
            "observer_free_scope_tag_v1",
        ],
        "observer_role": "downstream_only",
        "actual_nonzero_flow_discharged": f143["actual_nonzero_flow_discharged"],
        "barrier_protected_sign_discharged": f141["barrier_protected_sign_discharged"],
        "actual_observer_free_scope_discharged": f142["actual_observer_free_scope_discharged"],
        "actual_nontriviality_components_package_exported": actual_nontriviality_components_package_exported,
        "actual_nontriviality_assembly_witness_exported": False,
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
