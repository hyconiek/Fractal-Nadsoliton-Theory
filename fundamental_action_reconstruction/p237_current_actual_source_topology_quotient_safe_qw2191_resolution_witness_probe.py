#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F149 = ROOT / "generated" / "f149_first_actual_source_topology_quotient_safe_qw2191_resolution_witness_packet_summary.json"
OUT = ROOT / "generated" / "p237_current_actual_source_topology_quotient_safe_qw2191_resolution_witness_probe_summary.json"


def main() -> None:
    f149 = json.loads(IN_F149.read_text(encoding="utf-8"))

    passed = (
        f149["input_packet"] == "tau_src_candidate_v1"
        and f149["witness"] == "Phi_qw2191_safe_actual_witness_v1"
        and f149["codomain_target"] == "actual_quotient_safe_qw2191_resolution_target_v1"
        and f149["actual_basis_independent_selector_promotion_exported"] is True
        and f149["basis_independence_discharged"] is True
        and f149["kernel_alone_qw2191_obstruction_retained"] is True
        and f149["distinguished_quotient_class_exported"] is True
        and f149["raw_theta_uniqueness_claimed"] is False
        and f149["actual_quotient_safe_qw2191_resolution_exported"] is True
        and f149["qw2191_quotient_safe_discharged"] is True
        and f149["current_selector_closure"] is False
        and f149["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P237",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_WITNESS_BELOW_CURRENT_SELECTOR_CLOSURE_AFTER_P237"
            if passed
            else "P237_FAIL"
        ),
        "input_packet": f149["input_packet"],
        "witness": f149["witness"],
        "codomain_target": f149["codomain_target"],
        "support_packet_id": f149["support_packet_id"],
        "observer_role": f149["observer_role"],
        "tau_src_identified_with_s_prelm": f149["tau_src_identified_with_s_prelm"],
        "actual_basis_independent_selector_promotion_exported": f149["actual_basis_independent_selector_promotion_exported"],
        "basis_independence_discharged": f149["basis_independence_discharged"],
        "kernel_alone_qw2191_obstruction_retained": f149["kernel_alone_qw2191_obstruction_retained"],
        "distinguished_quotient_class_exported": f149["distinguished_quotient_class_exported"],
        "raw_theta_uniqueness_claimed": f149["raw_theta_uniqueness_claimed"],
        "actual_quotient_safe_qw2191_resolution_exported": f149["actual_quotient_safe_qw2191_resolution_exported"],
        "qw2191_quotient_safe_discharged": f149["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f149["current_selector_closure"],
        "current_global_qw2191_discharge": f149["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
