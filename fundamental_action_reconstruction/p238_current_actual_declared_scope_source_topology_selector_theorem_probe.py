#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F150 = ROOT / "generated" / "f150_first_actual_declared_scope_source_topology_selector_theorem_packet_summary.json"
OUT = ROOT / "generated" / "p238_current_actual_declared_scope_source_topology_selector_theorem_probe_summary.json"


def main() -> None:
    f150 = json.loads(IN_F150.read_text(encoding="utf-8"))

    passed = (
        f150["input_packet"] == "tau_src_candidate_v1"
        and f150["witness"] == "T14_src_selector_declared_scope_actual_witness_v1"
        and f150["codomain_target"] == "declared_scope_source_topology_selector_theorem_target_v1"
        and f150["actual_full_source_topology_nontriviality_witness_exported"] is True
        and f150["actual_observer_free_scope_discharged"] is True
        and f150["actual_selector_witness_exported"] is True
        and f150["basis_independence_discharged"] is True
        and f150["qw2191_quotient_safe_discharged"] is True
        and f150["tau_src_identified_with_s_prelm"] is False
        and f150["declared_scope_only"] is True
        and f150["admissible_strict_core_internal_selector_source_object_claimed"] is False
        and f150["declared_scope_source_topology_selector_theorem_exported"] is True
        and f150["current_selector_closure"] is False
        and f150["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P238",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_BELOW_CURRENT_SELECTOR_CLOSURE_AFTER_P238"
            if passed
            else "P238_FAIL"
        ),
        "input_packet": f150["input_packet"],
        "witness": f150["witness"],
        "codomain_target": f150["codomain_target"],
        "support_packet_id": f150["support_packet_id"],
        "observer_role": f150["observer_role"],
        "tau_src_identified_with_s_prelm": f150["tau_src_identified_with_s_prelm"],
        "declared_scope_only": f150["declared_scope_only"],
        "admissible_strict_core_internal_selector_source_object_claimed": f150["admissible_strict_core_internal_selector_source_object_claimed"],
        "actual_full_source_topology_nontriviality_witness_exported": f150["actual_full_source_topology_nontriviality_witness_exported"],
        "actual_observer_free_scope_discharged": f150["actual_observer_free_scope_discharged"],
        "actual_selector_witness_exported": f150["actual_selector_witness_exported"],
        "basis_independence_discharged": f150["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f150["qw2191_quotient_safe_discharged"],
        "declared_scope_source_topology_selector_theorem_exported": f150["declared_scope_source_topology_selector_theorem_exported"],
        "current_selector_closure": f150["current_selector_closure"],
        "current_global_qw2191_discharge": f150["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
