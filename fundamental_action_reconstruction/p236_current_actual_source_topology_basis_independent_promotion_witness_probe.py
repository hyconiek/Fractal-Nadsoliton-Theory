#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F148 = ROOT / "generated" / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
OUT = ROOT / "generated" / "p236_current_actual_source_topology_basis_independent_promotion_witness_probe_summary.json"


def main() -> None:
    f148 = json.loads(IN_F148.read_text(encoding="utf-8"))

    passed = (
        f148["input_packet"] == "tau_src_candidate_v1"
        and f148["witness"] == "Upsilon_sel_basis_actual_witness_v1"
        and f148["codomain_packet"] == "Sigma_sel_basis_free_target_v1"
        and f148["actual_basis_independent_selector_promotion_exported"] is True
        and f148["basis_free_axis_class_exported"] is True
        and f148["basis_free_signed_split_class_exported"] is True
        and f148["basis_free_scope_tag_exported"] is True
        and f148["observer_role"] == "downstream_only"
        and f148["tau_src_identified_with_s_prelm"] is False
        and f148["basis_independence_discharged"] is True
        and f148["qw2191_quotient_safe_discharged"] is False
        and f148["current_selector_closure"] is False
        and f148["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P236",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_BELOW_QW2191_AFTER_P236"
            if passed
            else "P236_FAIL"
        ),
        "input_packet": f148["input_packet"],
        "witness": f148["witness"],
        "codomain_packet": f148["codomain_packet"],
        "support_packet_id": f148["support_packet_id"],
        "observer_role": f148["observer_role"],
        "chart_bound_selector_input": f148["chart_bound_selector_input"],
        "tau_src_identified_with_s_prelm": f148["tau_src_identified_with_s_prelm"],
        "actual_full_source_topology_nontriviality_witness_exported": f148["actual_full_source_topology_nontriviality_witness_exported"],
        "full_source_topology_nontriviality_discharged": f148["full_source_topology_nontriviality_discharged"],
        "actual_selector_witness_exported": f148["actual_selector_witness_exported"],
        "basis_free_axis_class_exported": f148["basis_free_axis_class_exported"],
        "basis_free_signed_split_class_exported": f148["basis_free_signed_split_class_exported"],
        "basis_free_scope_tag_exported": f148["basis_free_scope_tag_exported"],
        "actual_basis_independent_selector_promotion_exported": f148["actual_basis_independent_selector_promotion_exported"],
        "basis_independence_discharged": f148["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f148["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f148["current_selector_closure"],
        "current_global_qw2191_discharge": f148["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
