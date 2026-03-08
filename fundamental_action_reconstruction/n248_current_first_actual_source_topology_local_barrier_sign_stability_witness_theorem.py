#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P228 = ROOT / "generated" / "p228_current_actual_source_topology_local_barrier_sign_stability_witness_probe_summary.json"
OUT = ROOT / "generated" / "n248_current_first_actual_source_topology_local_barrier_sign_stability_witness_theorem_summary.json"


def main() -> None:
    p228 = json.loads(IN_P228.read_text(encoding="utf-8"))

    passed = (
        p228["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_LOCAL_BARRIER_SIGN_STABILITY_WITNESS_BELOW_FULL_BARRIER_PROTECTED_SIGN_DISCHARGE_AFTER_P228"
        and p228["scalar_sign_component_from_f139"] == 1
        and p228["local_radius_positive"] is True
        and p228["interval_inside_positive_cos_domain"] is True
        and p228["barrier_protected_sign_discharged"] is False
        and p228["full_source_topology_nontriviality_discharged"] is False
        and p228["basis_independence_discharged"] is False
        and p228["qw2191_quotient_safe_discharged"] is False
        and p228["current_selector_closure"] is False
        and p228["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N248",
        "status": (
            "N248_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_LOCAL_BARRIER_SIGN_STABILITY_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N248_FAIL"
        ),
        "input_packet": p228["input_packet"],
        "witness": p228["witness"],
        "witness_definition": p228["witness_definition"],
        "local_radius": p228["local_radius"],
        "local_radius_definition": p228["local_radius_definition"],
        "local_radius_value": p228["local_radius_value"],
        "local_radius_positive": p228["local_radius_positive"],
        "scalar_sign_component_from_f139": p228["scalar_sign_component_from_f139"],
        "interval_inside_positive_cos_domain": p228["interval_inside_positive_cos_domain"],
        "local_interval_lower_endpoint": p228["local_interval_lower_endpoint"],
        "local_interval_upper_endpoint": p228["local_interval_upper_endpoint"],
        "observer_role": p228["observer_role"],
        "barrier_protected_sign_discharged": p228["barrier_protected_sign_discharged"],
        "full_source_topology_nontriviality_discharged": p228["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": p228["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p228["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p228["current_selector_closure"],
        "current_global_qw2191_discharge": p228["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
