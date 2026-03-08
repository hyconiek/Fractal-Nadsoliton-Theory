#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P227 = ROOT / "generated" / "p227_current_actual_source_topology_barrier_sign_component_witness_probe_summary.json"
OUT = ROOT / "generated" / "n247_current_first_actual_source_topology_barrier_sign_component_witness_theorem_summary.json"


def main() -> None:
    p227 = json.loads(IN_P227.read_text(encoding="utf-8"))

    passed = (
        p227["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_BELOW_FULL_BARRIER_PROTECTED_SIGN_DISCHARGE_AFTER_P227"
        and p227["nonzero_flow_component_support_from_f138"] is True
        and p227["barrier_margin_positive"] is True
        and p227["witness_value"] == 1
        and p227["barrier_protected_sign_discharged"] is False
        and p227["full_source_topology_nontriviality_discharged"] is False
        and p227["basis_independence_discharged"] is False
        and p227["qw2191_quotient_safe_discharged"] is False
        and p227["current_selector_closure"] is False
        and p227["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N247",
        "status": (
            "N247_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N247_FAIL"
        ),
        "input_packet": p227["input_packet"],
        "witness": p227["witness"],
        "witness_definition": p227["witness_definition"],
        "witness_value": p227["witness_value"],
        "barrier_margin": p227["barrier_margin"],
        "barrier_margin_definition": p227["barrier_margin_definition"],
        "barrier_margin_value": p227["barrier_margin_value"],
        "barrier_margin_positive": p227["barrier_margin_positive"],
        "nonzero_flow_component_support_from_f138": p227["nonzero_flow_component_support_from_f138"],
        "observer_role": p227["observer_role"],
        "barrier_protected_sign_discharged": p227["barrier_protected_sign_discharged"],
        "full_source_topology_nontriviality_discharged": p227["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p227["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p227["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p227["current_selector_closure"],
        "current_global_qw2191_discharge": p227["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
