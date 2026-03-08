#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P238 = ROOT / "generated" / "p238_current_actual_declared_scope_source_topology_selector_theorem_probe_summary.json"
OUT = ROOT / "generated" / "n258_current_first_declared_scope_source_topology_selector_theorem_summary.json"


def main() -> None:
    p238 = json.loads(IN_P238.read_text(encoding="utf-8"))

    passed = (
        p238["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_BELOW_CURRENT_SELECTOR_CLOSURE_AFTER_P238"
        and p238["declared_scope_only"] is True
        and p238["admissible_strict_core_internal_selector_source_object_claimed"] is False
        and p238["declared_scope_source_topology_selector_theorem_exported"] is True
        and p238["current_selector_closure"] is False
        and p238["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N258",
        "status": (
            "N258_DISCHARGED_CURRENT_FIRST_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_NO_FALSE_PASS"
            if passed
            else "N258_FAIL"
        ),
        "input_packet": p238["input_packet"],
        "witness": p238["witness"],
        "codomain_target": p238["codomain_target"],
        "support_packet_id": p238["support_packet_id"],
        "observer_role": p238["observer_role"],
        "tau_src_identified_with_s_prelm": p238["tau_src_identified_with_s_prelm"],
        "declared_scope_only": p238["declared_scope_only"],
        "admissible_strict_core_internal_selector_source_object_claimed": p238["admissible_strict_core_internal_selector_source_object_claimed"],
        "actual_full_source_topology_nontriviality_witness_exported": p238["actual_full_source_topology_nontriviality_witness_exported"],
        "actual_observer_free_scope_discharged": p238["actual_observer_free_scope_discharged"],
        "actual_selector_witness_exported": p238["actual_selector_witness_exported"],
        "basis_independence_discharged": p238["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p238["qw2191_quotient_safe_discharged"],
        "declared_scope_source_topology_selector_theorem_exported": p238["declared_scope_source_topology_selector_theorem_exported"],
        "current_selector_closure": p238["current_selector_closure"],
        "current_global_qw2191_discharge": p238["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
