#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F140 = ROOT / "generated" / "f140_first_actual_source_topology_local_barrier_sign_stability_witness_packet_summary.json"
OUT = ROOT / "generated" / "p228_current_actual_source_topology_local_barrier_sign_stability_witness_probe_summary.json"


def main() -> None:
    f140 = json.loads(IN_F140.read_text(encoding="utf-8"))

    passed = (
        f140["input_packet"] == "tau_src_candidate_v1"
        and f140["witness"] == "chi_src_local_barrier_sign_stability_witness_v1"
        and f140["scalar_sign_component_from_f139"] == 1
        and f140["local_radius_positive"] is True
        and f140["interval_inside_positive_cos_domain"] is True
        and f140["actual_local_sign_stability_witness_exported"] is True
        and f140["observer_role"] == "downstream_only"
        and f140["barrier_protected_sign_discharged"] is False
        and f140["full_source_topology_nontriviality_discharged"] is False
        and f140["basis_independence_discharged"] is False
        and f140["qw2191_quotient_safe_discharged"] is False
        and f140["current_selector_closure"] is False
        and f140["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P228",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_LOCAL_BARRIER_SIGN_STABILITY_WITNESS_BELOW_FULL_BARRIER_PROTECTED_SIGN_DISCHARGE_AFTER_P228"
            if passed
            else "P228_FAIL"
        ),
        "input_packet": f140["input_packet"],
        "witness": f140["witness"],
        "witness_definition": f140["witness_definition"],
        "local_radius": f140["local_radius"],
        "local_radius_definition": f140["local_radius_definition"],
        "local_radius_value": f140["local_radius_value"],
        "local_radius_positive": f140["local_radius_positive"],
        "scalar_sign_component_from_f139": f140["scalar_sign_component_from_f139"],
        "interval_inside_positive_cos_domain": f140["interval_inside_positive_cos_domain"],
        "local_interval_lower_endpoint": f140["local_interval_lower_endpoint"],
        "local_interval_upper_endpoint": f140["local_interval_upper_endpoint"],
        "observer_role": f140["observer_role"],
        "actual_local_sign_stability_witness_exported": f140["actual_local_sign_stability_witness_exported"],
        "barrier_protected_sign_discharged": f140["barrier_protected_sign_discharged"],
        "full_source_topology_nontriviality_discharged": f140["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f140["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f140["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f140["current_selector_closure"],
        "current_global_qw2191_discharge": f140["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
