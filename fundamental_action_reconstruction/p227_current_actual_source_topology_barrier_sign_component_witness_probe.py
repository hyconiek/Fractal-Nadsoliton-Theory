#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F139 = ROOT / "generated" / "f139_first_actual_source_topology_barrier_sign_component_witness_packet_summary.json"
OUT = ROOT / "generated" / "p227_current_actual_source_topology_barrier_sign_component_witness_probe_summary.json"


def main() -> None:
    f139 = json.loads(IN_F139.read_text(encoding="utf-8"))

    passed = (
        f139["input_packet"] == "tau_src_candidate_v1"
        and f139["witness"] == "psi_src_barrier_sign_component_witness_v1"
        and f139["nonzero_flow_component_support_from_f138"] is True
        and f139["barrier_margin_positive"] is True
        and f139["witness_value"] == 1
        and f139["actual_scalar_sign_component_witness_exported"] is True
        and f139["observer_role"] == "downstream_only"
        and f139["barrier_protected_sign_discharged"] is False
        and f139["full_source_topology_nontriviality_discharged"] is False
        and f139["basis_independence_discharged"] is False
        and f139["qw2191_quotient_safe_discharged"] is False
        and f139["current_selector_closure"] is False
        and f139["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P227",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_BELOW_FULL_BARRIER_PROTECTED_SIGN_DISCHARGE_AFTER_P227"
            if passed
            else "P227_FAIL"
        ),
        "input_packet": f139["input_packet"],
        "witness": f139["witness"],
        "witness_definition": f139["witness_definition"],
        "witness_value": f139["witness_value"],
        "barrier_margin": f139["barrier_margin"],
        "barrier_margin_definition": f139["barrier_margin_definition"],
        "barrier_margin_value": f139["barrier_margin_value"],
        "barrier_margin_positive": f139["barrier_margin_positive"],
        "nonzero_flow_component_support_from_f138": f139["nonzero_flow_component_support_from_f138"],
        "observer_role": f139["observer_role"],
        "actual_scalar_sign_component_witness_exported": f139["actual_scalar_sign_component_witness_exported"],
        "barrier_protected_sign_discharged": f139["barrier_protected_sign_discharged"],
        "full_source_topology_nontriviality_discharged": f139["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f139["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f139["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f139["current_selector_closure"],
        "current_global_qw2191_discharge": f139["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
