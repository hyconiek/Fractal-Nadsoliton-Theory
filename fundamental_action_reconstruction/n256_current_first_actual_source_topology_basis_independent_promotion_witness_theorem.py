#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P236 = ROOT / "generated" / "p236_current_actual_source_topology_basis_independent_promotion_witness_probe_summary.json"
OUT = ROOT / "generated" / "n256_current_first_actual_source_topology_basis_independent_promotion_witness_theorem_summary.json"


def main() -> None:
    p236 = json.loads(IN_P236.read_text(encoding="utf-8"))

    passed = (
        p236["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_BELOW_QW2191_AFTER_P236"
        and p236["actual_basis_independent_selector_promotion_exported"] is True
        and p236["basis_independence_discharged"] is True
        and p236["qw2191_quotient_safe_discharged"] is False
        and p236["current_selector_closure"] is False
        and p236["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N256",
        "status": (
            "N256_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N256_FAIL"
        ),
        "input_packet": p236["input_packet"],
        "witness": p236["witness"],
        "codomain_packet": p236["codomain_packet"],
        "support_packet_id": p236["support_packet_id"],
        "observer_role": p236["observer_role"],
        "chart_bound_selector_input": p236["chart_bound_selector_input"],
        "tau_src_identified_with_s_prelm": p236["tau_src_identified_with_s_prelm"],
        "actual_full_source_topology_nontriviality_witness_exported": p236["actual_full_source_topology_nontriviality_witness_exported"],
        "full_source_topology_nontriviality_discharged": p236["full_source_topology_nontriviality_discharged"],
        "actual_selector_witness_exported": p236["actual_selector_witness_exported"],
        "basis_free_axis_class_exported": p236["basis_free_axis_class_exported"],
        "basis_free_signed_split_class_exported": p236["basis_free_signed_split_class_exported"],
        "basis_free_scope_tag_exported": p236["basis_free_scope_tag_exported"],
        "actual_basis_independent_selector_promotion_exported": p236["actual_basis_independent_selector_promotion_exported"],
        "basis_independence_discharged": p236["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p236["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p236["current_selector_closure"],
        "current_global_qw2191_discharge": p236["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
