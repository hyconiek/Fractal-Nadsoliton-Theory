#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P235 = ROOT / "generated" / "p235_current_actual_source_topology_selector_witness_probe_summary.json"
OUT = ROOT / "generated" / "n255_current_first_actual_source_topology_selector_witness_theorem_summary.json"


def main() -> None:
    p235 = json.loads(IN_P235.read_text(encoding="utf-8"))

    passed = (
        p235["status"]
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_BELOW_BASIS_INDEPENDENCE_AND_QW2191_AFTER_P235"
        and p235["actual_selector_witness_exported"] is True
        and p235["chart_bound_preobserver_realization"] is True
        and p235["tau_src_identified_with_s_prelm"] is False
        and p235["basis_independence_discharged"] is False
        and p235["qw2191_quotient_safe_discharged"] is False
        and p235["current_selector_closure"] is False
        and p235["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N255",
        "status": (
            "N255_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_THEOREM_NO_FALSE_PASS"
            if passed
            else "N255_FAIL"
        ),
        "input_packet": p235["input_packet"],
        "witness": p235["witness"],
        "codomain_packet": p235["codomain_packet"],
        "support_packet_id": p235["support_packet_id"],
        "observer_role": p235["observer_role"],
        "chart_bound_preobserver_realization": p235["chart_bound_preobserver_realization"],
        "tau_src_identified_with_s_prelm": p235["tau_src_identified_with_s_prelm"],
        "actual_full_source_topology_nontriviality_witness_exported": p235["actual_full_source_topology_nontriviality_witness_exported"],
        "full_source_topology_nontriviality_discharged": p235["full_source_topology_nontriviality_discharged"],
        "actual_selector_witness_exported": p235["actual_selector_witness_exported"],
        "basis_independence_discharged": p235["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p235["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p235["current_selector_closure"],
        "current_global_qw2191_discharge": p235["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
