#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P237 = ROOT / "generated" / "p237_current_actual_source_topology_quotient_safe_qw2191_resolution_witness_probe_summary.json"
OUT = ROOT / "generated" / "n257_current_first_actual_source_topology_quotient_safe_qw2191_resolution_witness_theorem_summary.json"


def main() -> None:
    p237 = json.loads(IN_P237.read_text(encoding="utf-8"))

    passed = (
        p237["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_WITNESS_BELOW_CURRENT_SELECTOR_CLOSURE_AFTER_P237"
        and p237["actual_quotient_safe_qw2191_resolution_exported"] is True
        and p237["qw2191_quotient_safe_discharged"] is True
        and p237["current_selector_closure"] is False
        and p237["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N257",
        "status": (
            "N257_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N257_FAIL"
        ),
        "input_packet": p237["input_packet"],
        "witness": p237["witness"],
        "codomain_target": p237["codomain_target"],
        "support_packet_id": p237["support_packet_id"],
        "observer_role": p237["observer_role"],
        "tau_src_identified_with_s_prelm": p237["tau_src_identified_with_s_prelm"],
        "actual_basis_independent_selector_promotion_exported": p237["actual_basis_independent_selector_promotion_exported"],
        "basis_independence_discharged": p237["basis_independence_discharged"],
        "kernel_alone_qw2191_obstruction_retained": p237["kernel_alone_qw2191_obstruction_retained"],
        "distinguished_quotient_class_exported": p237["distinguished_quotient_class_exported"],
        "raw_theta_uniqueness_claimed": p237["raw_theta_uniqueness_claimed"],
        "actual_quotient_safe_qw2191_resolution_exported": p237["actual_quotient_safe_qw2191_resolution_exported"],
        "qw2191_quotient_safe_discharged": p237["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p237["current_selector_closure"],
        "current_global_qw2191_discharge": p237["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
